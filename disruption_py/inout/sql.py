#!/usr/bin/env python3

"""
Module for managing SQL database connections.
"""

import os
import threading
from typing import List
from urllib.parse import quote_plus

import numpy as np
import pandas as pd
import pyodbc
from loguru import logger
from sqlalchemy import create_engine

from disruption_py.config import config
from disruption_py.core.utils.misc import shot_log_msg, without_duplicates
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
        conn_str = ";".join([f"{k}={v}" for k, v in params.items()])
        return conn_str

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
                "Connecting to database for thread {current_thread}",
                current_thread=str(current_thread),
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
            logger.error("Query failed with error {e}, returning None", e=e)
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

    def add_shot_data(
        self,
        shot_id: int,
        shot_data: pd.DataFrame,
        update=False,
        override_columns: List[str] = None,
    ):
        """
        Upload shot to SQL database.

        Either inserts or updates shot data depending on whether a shot already exists
        in database. If shot exists, then the timebase of the shot data must match
        the timebase of the shot in the database.

        Parameters
        ----------
        shot_id : int
            Shot id of the shot being modified
        shot_data : pd.DataFrame
            Dataframe containing shot data for update
        update : bool
            Whether to update shot data if the shot already exists in database.
            Update will happen regardless of whether the column being updated is
            all nil. Default value is False.
        override_columns : List[str]
            List of protected columns that can still be updated. Update must be
            true for input values in the columns to be changed. Default value is [].
        """
        if self.write_database_table_name is None:
            raise ValueError(
                "specify write_database_table_name in the configuration before "
                + "adding shot data"
            )
        curr_df = pd.read_sql_query(
            f"select * from {self.write_database_table_name} where shot={shot_id} "
            + "order by time",
            self.engine,
        )

        if len(curr_df) == 0:
            return self._insert_shot_data(
                curr_df=curr_df,
                shot_data=shot_data,
                table_name=self.write_database_table_name,
            )
        if (
            len(curr_df) == len(shot_data)
            and (
                (curr_df["time"] - shot_data["time"]).abs() < config().time_const
            ).all()
        ):
            return self._update_shot_data(
                shot_id=shot_id,
                curr_df=curr_df,
                shot_data=shot_data,
                update=update,
                table_name=self.write_database_table_name,
                override_columns=override_columns,
            )

        logger.error("Invalid timebase for data output")
        return False

    def _insert_shot_data(
        self,
        curr_df: pd.DataFrame,
        shot_data: pd.DataFrame,
        table_name: str,
    ):
        """
        Insert shot data into SQL table.

        Assumes that the shot id does not already exist in the database.
        """

        identity_column_names = self._get_identity_column_names(table_name)

        matching_columns_shot_data = pd.DataFrame()
        for column_name in curr_df.columns:
            if column_name in identity_column_names:
                continue

            if column_name in shot_data.columns:
                matching_columns_shot_data[column_name] = shot_data[column_name]

        matching_columns_shot_data = matching_columns_shot_data.replace({np.nan: None})

        column_names = matching_columns_shot_data.columns.tolist()
        sql_column_names = ", ".join(column_names)
        parameter_markers = "(" + ", ".join(["?"] * len(column_names)) + ")"
        with self.conn.cursor() as curs:
            data_tuples = list(
                matching_columns_shot_data.itertuples(index=False, name=None)
            )
            curs.executemany(
                f"insert into {table_name} ({sql_column_names}) values "
                + f"{parameter_markers}",
                data_tuples,
            )
        return True

    def _update_shot_data(
        self,
        shot_id: int,
        curr_df: pd.DataFrame,
        shot_data: pd.DataFrame,
        update: bool,
        table_name: str,
        override_columns: List[str] = None,
    ):
        """
        Update shot data into SQL table.

        Assumes that the shot id already exist in the database and the timebase of
        shot_data is the same as curr_df.

        Parameters
        ----------
        curr_df : pd.DataFrame
            Data currently in sql database.
        shot_data : pd.DataFrame
            Dataframe containing shot data for update.
        update : bool
            Whether to update shot data if the shot already exists in database.
            Update will happen regardless of whether the column being updated is
            all nil. Default value is False.
        override_columns : List[str]
            List of columns that can should still be updated. Update must be true
            for input values in the columns to be changed. Default value is [].
        table_name : str
            Name of the table for data insert or update. Default value is
            "disruption_warning".
        """
        override_columns = override_columns or []

        update_columns_shot_data = pd.DataFrame()
        for column_name in curr_df.columns:
            if column_name in config().inout.sql.protected_columns or (
                column_name in self.protected_columns
                and column_name not in override_columns
            ):
                continue

            if (
                column_name in shot_data.columns
                and not shot_data[column_name].isna().all()
                and (update or curr_df[column_name].isna().all())
            ):
                update_columns_shot_data[column_name] = shot_data[column_name]
        # pyodbc will fill SQL with NULL for None, but not for np.nan
        update_columns_shot_data = update_columns_shot_data.replace({np.nan: None})
        ko_rows = 0
        with self.conn.cursor() as curs:
            for index, row in enumerate(
                update_columns_shot_data.itertuples(index=False, name=None)
            ):
                update_column_names = list(update_columns_shot_data.columns)
                sql_set_string = ", ".join(
                    [f"{col} = ?" for col in update_column_names]
                )
                sql_command = (
                    f"UPDATE {table_name} SET {sql_set_string} "
                    + "WHERE time BETWEEN ? AND ? AND shot = ?;"
                )
                t = curr_df["time"][index] + np.array([-0.5, 0.5]) * config().time_const
                curs.execute(sql_command, row + (*t, str(shot_id)))
                ko_rows += curs.rowcount == 0
        if ko_rows:
            logger.error(shot_log_msg(shot_id, f"Could not update {ko_rows} rows."))
        return True

    def _get_identity_column_names(self, table_name: str):
        """Get which column names are identity columns in table."""
        queries = {
            "mssql": """
                     SELECT
                       c.name AS ColumnName
                     FROM
                       sys.columns c
                       INNER JOIN sys.tables t ON c.object_id = t.object_id
                       LEFT JOIN sys.identity_columns ic ON ic.object_id = c.object_id
                       AND ic.column_id = c.column_id
                     WHERE
                       t.name = '{table_name}'
                       AND ic.object_id IS NOT NULL
                     """,
            "mysql": """
                     SELECT
                       COLUMN_NAME AS ColumnName
                     FROM
                       INFORMATION_SCHEMA.COLUMNS
                     WHERE
                       TABLE_NAME = '{table_name}'
                       AND EXTRA LIKE '%auto_increment%';
                     """,
        }
        with self.conn.cursor() as curs:
            query = queries[self.dialect].format(table_name=table_name)
            logger.trace("Executing query: {query}", query=query)
            curs.execute(query)
            return [row[0] for row in curs.fetchall()]

    def remove_shot_data(self, shotlist: List[int]):
        """Remove shot data from the test SQL table."""
        table_name = self.write_database_table_name
        if not table_name.endswith("_test"):
            raise ValueError("Deletion is restricted to tables ending in '_test'.")
        shots = [str(s) for s in shotlist]
        with self.conn.cursor() as curs:
            query = f"delete from {table_name} where shot in ({', '.join(shots)})"
            logger.debug("Executing query: '{query}'", query=query)
            curs.execute(query)
            logger.debug("Deleted: {rows} rows", rows=curs.rowcount)
            return curs.rowcount > 0

    def add_column(self, col_name, var_type="TEXT"):
        """Add column to SQL table without filling in data for column."""
        if self.write_database_table_name is None:
            raise ValueError(
                "specify write_database_table_name in the configuration before "
                + "adding shot data"
            )
        self.query(
            f"alter table {self.write_database_table_name} add {col_name} {var_type};",
            use_pandas=False,
        )
        return True

    def remove_column(self, col_name):
        """Remove column from SQL table"""
        if self.write_database_table_name is None:
            raise ValueError(
                "specify write_database_table_name in the configuration before "
                + "adding shot data"
            )
        if col_name in self.protected_columns:
            logger.error(
                "Failed to drop protected column {col_name}", col_name=col_name
            )
            return False
        self.query(
            f"alter table {self.write_database_table_name} drop column {col_name};",
            use_pandas=False,
        )
        return True

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

    def get_disruption_shotlist(self):
        """
        Get Pandas DataFrame of all disruptive shots and times from the disruption
        table. Can be set as a cross-reference to determine whether a given shot
        is disruptive or not (all shots in this table are disruptive) and contain
        a t_disrupt.
        """
        return self.query("select distinct shot from disruptions order by shot")

    def get_disruption_warning_shotlist(self):
        """
        Get Pandas DataFrame of all shots in the disruption_warning table. NOTE:
        The disruption_warning table contains ONLY a subset of shots in this table
        """
        return self.query("select distinct shot from disruption_warning order by shot")


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

    def get_disruption_shotlist(self, **_kwargs):
        return []

    def get_disruption_warning_shotlist(self, **_kwargs):
        return []
