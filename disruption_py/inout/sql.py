#!/usr/bin/env python3

"""
Module for managing SQL database connections.
"""

import os
import threading
from typing import List
from urllib.parse import quote_plus

import pandas as pd
import pyodbc
from loguru import logger
from sqlalchemy import create_engine

from disruption_py.config import config
from disruption_py.core.utils.misc import without_duplicates
from disruption_py.core.utils.shared_instance import SharedInstance
from disruption_py.machine.tokamak import Tokamak


class ShotDatabase:
    """
    Handles grabbing data from MySQL server.
    """

    def __init__(
        self,
        driver,
        host,
        port,
        db_name,
        user,
        passwd,
        protected_columns=None,
        write_database_table_name=None,
        **_kwargs,
    ):

        if protected_columns is None:
            protected_columns = []

        logger.debug(
            "Database initialization: {user}@{host}/{db_name}",
            user=user,
            host=host,
            db_name=db_name,
        )
        drivers = pyodbc.drivers()
        if driver in drivers:
            # exact driver
            self.driver = driver
        elif any(d.startswith(driver) for d in drivers):
            # fuzzy driver
            self.driver = next(d for d in drivers if d.startswith(driver))
            logger.info(
                "Database driver fallback: '{driver}' -> '{class_driver}'",
                driver=driver,
                class_driver=self.driver,
            )
        else:
            # first driver
            self.driver = drivers[0]
            logger.warning(
                "Database driver fallback: '{driver}' -> '{class_driver}'",
                driver=driver,
                class_driver=self.driver,
            )
        self.host = host
        self.port = port
        self.db_name = db_name
        self.user = user
        self.passwd = passwd
        self.protected_columns = protected_columns
        self.write_database_table_name = write_database_table_name

        self.dialect = "mysql" if "mysql" in self.driver.lower() else "mssql"
        self.connection_string = self._get_connection_string(self.db_name)
        self._thread_connections = {}
        quoted_connection_string = quote_plus(self.connection_string)
        self.engine = create_engine(
            f"{self.dialect}+pyodbc:///?odbc_connect={quoted_connection_string}"
        )

    @classmethod
    def from_config(cls, tokamak: Tokamak):
        """
        Initialize database from config.
        """

        db_conf = config(tokamak).inout.sql

        # read sybase login
        if any(f"db_{key}" not in db_conf for key in ["user", "pass"]):
            db_name = db_conf["db_name"]
            for name in [db_name.lower(), db_name.upper()]:
                profile = os.path.expanduser(f"~/{name}.sybase_login")
                if not os.path.exists(profile):
                    continue
                with open(profile, "r", encoding="utf-8") as fio:
                    db_conf["db_user"], db_conf["db_pass"] = fio.read().split()[-2:]
                break
            else:
                raise ValueError("could not read DB username and password.")

        return SharedInstance(ShotDatabase).get_instance(
            driver=db_conf["driver"],
            host=db_conf["host"],
            port=db_conf["port"],
            db_name=db_conf["db_name"],
            user=db_conf["db_user"],
            passwd=db_conf["db_pass"],
            protected_columns=without_duplicates(db_conf["protected_columns"]),
            write_database_table_name=db_conf.get("write_database_table_name"),
        )

    def _get_connection_string(self, db_name):
        params = {
            "DRIVER": self.driver,
            "SERVER": self.host,
            "PORT": self.port,
            "DATABASE": db_name,
            "UID": self.user,
            "PWD": self.passwd,
            "TrustServerCertificate": "yes",
            "Connection Timeout": 60,
        }
        if self.driver.lower().startswith("odbc"):
            params["SERVER"] += f",{params.pop('PORT')}"
        return ";".join([f"{k}={v}" for k, v in params.items()])

    @property
    def conn(self):
        """
        Property returning a connection to sql database.

        If a connection exists for the given thread returns that connection,
        otherwise creates a new connection

        Returns
        -------
        _type_
            Database connection
        """
        current_thread = threading.current_thread()
        if current_thread not in self._thread_connections:
            logger.debug(
                "PID #{pid} | Connecting to SQL database: -> {server}",
                pid=threading.get_native_id(),
                server=self.host,
            )
            self._thread_connections[current_thread] = pyodbc.connect(
                self.connection_string
            )
        return self._thread_connections[current_thread]

    def query(self, query: str, use_pandas=True):
        """
        query sql database

        Parameters
        ----------
        query : str
            The query string
        use_pandas : bool, optional
            Whether pd.read_sql_query should be used to run the query. Default value
            is true.

        Returns
        -------
        Any
            Result of query
        """
        if "alter" in query.lower():
            if query.lower() in self.protected_columns:
                return 0
        elif use_pandas:
            return pd.read_sql_query(query, self.engine)
        curs = self.conn.cursor()
        output = None
        try:
            curs.execute(query)
            if "select" in query.lower():
                output = curs.fetchall()
        except pyodbc.DatabaseError as e:
            logger.error("Query failed: {e}", e=repr(e))
            logger.opt(exception=True).debug(e)
        curs.close()
        return output

    def get_version(self):
        """
        Query the version of the SQL database.
        """
        mysql = "mysql" in self.driver.lower()
        query = "select " + ("version()" if mysql else "@@version")
        version = self.query(query, use_pandas=False)
        return version[0][0]

    def get_shots_data(
        self,
        shotlist: List[int],
        cols: List[str] = None,
        sql_table="disruption_warning",
    ):
        """
        get_shots_data retrieves columns from sql data for given shotlist

        Parameters
        ----------
        shotlist : List[int]
            List of shot ids to get data for.
        cols : List[str]
            List of columns to retrieve. Default value is ["*"], meaning all columns.
        sql_table : str, optional
            The sql_table to retrieve data from. Default value is "disruption_warning".

        Returns
        -------
        pd.Dataframe
            Dataframe containing queried data
        """
        if cols is None:
            cols = ["*"]
        cols = ", ".join(str(col) for col in cols)
        shotlist = ",".join(str(shot) for shot in shotlist)
        query = f"select {cols} from {sql_table}"
        if shotlist is None:
            query += " order by time"
        else:
            query += f" where shot in ({shotlist}) order by shot, time"
        shot_df = pd.read_sql_query(query, self.engine)
        shot_df.columns = shot_df.columns.str.lower()
        return shot_df

    def get_disruption_time(self, shot_id):
        """
        Get disruption time for shot_id or None if there was no disruption.
        """
        with self.conn.cursor() as curs:
            curs.execute(f"select t_disrupt from disruptions where shot = {shot_id}")
            t_disrupt = curs.fetchall()
        if len(t_disrupt) == 0:
            return None
        t_disrupt = t_disrupt[0][0]
        return t_disrupt


class DummyObject:
    """
    A dummy connection object.
    """

    def __getattr__(self, name):
        # Return self for any attribute or method call
        return self

    def __call__(self, *args, **kwargs):
        # Return self for any method call
        return self


class DummyDatabase(ShotDatabase):
    """
    A database class that does not require connecting to an SQL server and returns
    no data.
    """

    # pylint: disable-next=super-init-not-called
    def __init__(self, **kwargs):
        pass

    @classmethod
    # pylint: disable-next=missing-function-docstring
    def initializer(cls, **_kwargs):
        return cls()

    @property
    def conn(self):
        return DummyObject()

    # pylint: disable-next=arguments-differ
    def query(self, **_kwargs):
        return pd.DataFrame()

    # pylint: disable-next=arguments-differ
    def get_shots_data(self, **_kwargs):
        return pd.DataFrame()

    # pylint: disable-next=arguments-differ
    def get_disruption_time(self, **_kwargs):
        return None
