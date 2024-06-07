#!/usr/bin/env python3

from dataclasses import Field, dataclass
import logging
import os
import threading
from typing import List
from urllib.parse import quote_plus

import numpy as np
import pandas as pd
import pyodbc
from sqlalchemy import Engine, create_engine

from disruption_py.utils.constants import BASE_PROTECTED_COLUMNS, TIME_CONST, DATABASE_CONSTANTS
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment

@dataclass
class DBConnection:
    """
    Class for holding database connection information.
    """
    logger = logging.getLogger("disruption_py")
    connection_string : str
    engine : Engine = Field(init=False)
    _thread_connections : dict = {}
    
    def __post_init__(self):
        quoted_connection_string = quote_plus(self.connection_string)
        self.engine = create_engine(f"mssql+pyodbc:///?odbc_connect={quoted_connection_string}")
        
    @property
    def conn(self):
        """Property returning a connection to sql database.

        If a connection exists for the given thread returns that connection, otherwise creates a new connection

        Returns
        -------
        _type_
            Database connection
        """
        current_thread = threading.current_thread()
        if current_thread not in self._thread_connections:
            self.logger.info(f"Connecting to database for thread {current_thread}")
            self._thread_connections[current_thread] = pyodbc.connect(
                self.connection_string
            )
        return self._thread_connections[current_thread]

class ShotDatabase:
    """
    Handles grabbing data from MySQL server.
    """

    logger = logging.getLogger("disruption_py")

    def __init__(
        self, driver, host, port, db_names, user, passwd, default_db=None, protected_columns=[], **kwargs
    ):
        self.logger.info(f"Database initialization: {user}@{host}/{db_names}")
        drivers = pyodbc.drivers()
        if driver in drivers:
            self.driver = driver
        else:
            self.driver = drivers[0]
            self.logger.warning(
                f"Database driver fallback: '{driver}' -> '{self.driver}'"
            )
        self.host = host
        self.port = port
        self.user = user
        self.passwd = passwd
        self.protected_columns = protected_columns
        
        if default_db is None:
            default_db = db_names[0]
        if default_db not in db_names:
            raise ValueError(f"Default database {default_db} not in db_names: {db_names}")
        self.default_db = default_db
        
        self.dbs : dict[str, DBConnection] = {}
        for db_name in db_names:
            connection_string = self._get_connection_string(db_name)
            self.dbs[db_name] = DBConnection(connection_string)
            
    def from_config(self, tokamak : Tokamak=None):
        """
        Initialize database from config file.
        """
        if tokamak is None:
            tokamak = get_tokamak_from_environment()
            
        constants = DATABASE_CONSTANTS[tokamak.value]
        profile_path = constants["profile_path"]
        profile = os.path.expanduser(profile_path)
        with open(profile, "r") as fio:
            db_user, db_pass = fio.read().split()[-2:]
        
        return ShotDatabase(
            **constants,
            user=db_user,
            passwd=db_pass,
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
        if "ODBC" in self.driver:
            params["SERVER"] += f",{params.pop('PORT')}"
        conn_str = ";".join([f"{k}={v}" for k, v in params.items()])
        return conn_str

    def query(self, query: str, use_pandas=True, db_name=None):
        """query sql database

        Parameters
        ----------
        query : str
            The query string
        use_pandas : bool, optional
            Whether pd.read_sql_query should be used to run the query. Default value is true.

        Returns
        -------
        Any
            Result of query
        """
        if db_name is None:
            db_name = self.default_db
        
        if "alter" in query.lower():
            if query.lower() in self.protected_columns:
                return 0
        elif use_pandas and "select" in query.lower():
            return pd.read_sql_query(query, self.dbs[db_name].engine)
        curs = self.dbs[db_name].conn.cursor()
        output = None
        try:
            curs.execute(query)
            if "select" in query.lower():
                output = curs.fetchall()
        except pyodbc.DatabaseError as e:
            print(e)
            self.logger.debug(e)
            self.logger.error("Query failed, returning None")
        curs.close()
        return output

    def add_shot_data(
        self,
        shot_id: int,
        shot_data: pd.DataFrame,
        update=False,
        override_columns: List[str] = None,
        db_name=None,
        table_name="disruption_warning",
    ):
        """
        Upload shot to SQL database.

        Either inserts or updates shot data depending on whether shot already exists in database.
        If shot exists, then thye timebase of the shot data must match the timebase of the shot in the database.

        Parameters
        ----------
        shot_id : int
            Shot id of the shot being modified
        shot_data : pd.DataFrame
            Dataframe containing shot data for update
        update : bool
            Whether to update shot data if the shot already exists in database. Update will happen regardless if
            the column being updated is all nil. Default value is False.
        override_columns : List[str]
            List of protecrted columns that can should still be updated. Update must be true for existing values in the
            columns to be changed. Default value is [].
        table_name : str
            Name of the table for data insert or update. Default value is "disruption_warning".
        """
        if db_name is None:
            db_name = self.default_db
            
        curr_df = pd.read_sql_query(
            f"select * from {table_name} where shot in {shot_id} order by time",
            self.dbs[db_name].engine,
        )

        if len(curr_df == 0):
            return self._insert_shot_data(curr_df, shot_data)
        elif (
            len(curr_df) == len(shot_data)
            and ((curr_df["time"] - shot_data["time"]).abs() < TIME_CONST).all()
        ):
            return self._update_shot_data(curr_df, shot_data, update, override_columns)

        self.logger.error("Invalid timebase for data output")
        return False

    def _insert_shot_data(
        self,
        curr_df: pd.DataFrame,
        shot_data: pd.DataFrame,
        db_name=None,
        table_name="disruption_warning",
    ):
        """
        Insert shot data into SQL table.

        Assumes that the shot id does not already exist in the database.
        """
        if db_name is None:
            db_name = self.default_db
        matching_columns_shot_data = pd.DataFrame()
        for column_name in curr_df.columns:
            if column_name in shot_data.columns:
                matching_columns_shot_data[column_name] = shot_data[column_name]
            else:
                matching_columns_shot_data[column_name] = np.nan
        with self.dbs[db_name].conn.cursor() as curs:
            data_tuples = list(
                matching_columns_shot_data.itertuples(index=False, name=None)
            )
            curs.executemany(f"""insert into {table_name} values""", data_tuples)
        return True

    def _update_shot_data(
        self,
        curr_df: pd.DataFrame,
        shot_data: pd.DataFrame,
        update: bool,
        override_columns: List[str] = None,
        db_name=None,
        table_name="disruption_warning",
    ):
        """
        Update shot data into SQL table.

        Assumes that the shot id already exist in the database and timebase of shot_data is the same as curr_df.

        Parameters
        ----------
        curr_df : pd.DataFrame
            Data currently in sql database.
        shot_data : pd.DataFrame
            Dataframe containing shot data for update.
        update : bool
            Whether to update shot data if the shot already exists in database. Update will happen regardless if
            the column being updated is all nil. Default value is False.
        override_columns : List[str]
            List of protecrted columns that can should still be updated. Update must be true for existing values in the
            columns to be changed. Default value is [].
        table_name : str
            Name of the table for data insert or update. Default value is "disruption_warning".
        """
        if db_name is None:
            db_name = self.default_db
        override_columns = override_columns or []

        update_columns_shot_data = pd.DataFrame()
        for column_name in curr_df.columns:
            if column_name in BASE_PROTECTED_COLUMNS or (
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
        with self.dbs[db_name].conn.cursor() as curs:
            for index, row in update_columns_shot_data.iterrows():
                update_columns = update_columns_shot_data.columns
                sql_command = f"UPDATE {table_name} SET {', '.join([f'{col} = ?' for col in update_columns])} WHERE time = ?;"
                curs.execute(sql_command, row + curr_df["time"][index])
        return True

    def remove_shot_data(self, shot_id, db_name=None, table_name="disruption_warning"):
        """Remove shot from SQL table."""
        if db_name is None:
            db_name = self.default_db
        data_df = pd.read_sql_query(
            f"""select * from {table_name} where shot = {shot_id} order by time""",
            self.dbs[db_name].engine,
        )
        if len(data_df) == 0:
            self.logger.info(f"Shot {shot_id} does not exist in database")
            return False
        with self.dbs[db_name].conn.cursor() as curs:
            curs.execute(f"delete from disruption_warning where shot = {shot_id}")
        return True

    def add_column(self, col_name, var_type="TEXT", db_name=None, table_name="disruption_warning"):
        """Add column to SQL table without filling in data for column."""
        try:
            self.query(
                f"alter table {table_name} add {col_name} {var_type};", use_pandas=False, db_name=db_name,
            )
            return True
        except Exception as e:
            self.logger.error(f"Failed to add column {col_name} with error {e}")
            return False

    def remove_column(self, col_name, db_name=None, table_name="disruption_warning"):
        """Remove column from SQL table"""
        if col_name in self.protected_columns:
            self.logger.error(f"Failed to drop protected column {col_name}")
            return False
        try:
            self.query(
                f"alter table {table_name} drop column {col_name};", use_pandas=False, db_name=db_name,
            )
            return True
        except Exception as e:
            self.logger.error(f"Failed to drop column {col_name} with error {e}")
            return False

    def get_shots_data(
        self,
        shot_ids: List[int],
        cols: List[str] = ["*"],
        db_name=None,
        sql_table="disruption_warning",
    ):
        """get_shots_data retrieves columns from sql data for given shot_ids

        Parameters
        ----------
        shot_ids : List[int]
            List of shot ids to get data for.
        cols : List[str]
            List of columns to retrieve. Default value is ["*"], meaning all columns.
        sql_table : str, optional
            The sql_table to retrieve data from. Default value is "disruption_warning".

        Returns
        -------
        pd.Dataframe
            Dataframe containing querried data
        """
        if db_name is None:
            db_name = self.default_db
        shot_ids = ",".join([str(shot_id) for shot_id in shot_ids])
        selected_cols = f"{cols[0]}"
        if len(cols) > 1:
            selected_cols += "".join([f", {col}" for col in cols[1:]])
        if shot_ids is None:
            query = f"select {selected_cols} from {sql_table} order by time"
        else:
            query = f"select {selected_cols} from {sql_table} where shot in ({shot_ids}) order by time"
        shot_df = pd.read_sql_query(query, self.dbs[db_name].engine)
        return shot_df

    def get_disruption_time(self, shot_id, db_name=None,):
        """
        Get disruption time for shot_id or None if there was no disruption.
        """
        if db_name is None:
            db_name = self.default_db
        with self.dbs[db_name].conn.cursor() as curs:
            curs.execute(f"select t_disrupt from disruptions where shot = {shot_id}")
            t_disrupt = curs.fetchall()
        if len(t_disrupt) == 0:
            return None
        t_disrupt = t_disrupt[0][0]
        return t_disrupt

    def get_disruption_shotlist(self, db_name=None):
        """
        Get pandas dataframe of all disruptive shots and times from the disruption table. Can be sed as a cross-reference to determine whether a given shot is disruptive or not (all shots in this table are disruptive) and contain a t_disrupt.
        """
        return self.query("select distinct shot from disruptions order by shot", db_name=db_name)

    def get_disruption_warning_shotlist(self, db_name=None):
        """
        Get pandas dataframe of all shots in the disruption_warning table. NOTE: The disruption_warning table contains ONLY a subset of shots in this table
        """
        return self.query("select distinct shot from disruption_warning order by shot", db_name=db_name)


class DummyObject:
    def __getattr__(self, name):
        # Return self for any attribute or method call
        return self

    def __call__(self, *args, **kwargs):
        # Return self for any method call
        return self


class DummyDatabase(ShotDatabase):
    """
    A database class that does not require connecting to an SQL server but returns no data.

    Note: On CMod, disruption time data and any derrivative values will not be correct

    Examples
    --------
    >>> cmod_handler = CModHandler(database_initializer=DummyDatabase.default)
    >>> shot_data = cmod_handler.get_shots_data(shot_ids_request=[1150805012])
    <pd.DataFrame>
    """

    def __init__(self, **kwargs):
        pass

    @classmethod
    def default(cls, **kwargs):
        return cls()

    @property
    def conn(self, **kwargs):
        return DummyObject()

    def query(self, **kwargs):
        return pd.DataFrame()

    def get_shots_data(sefl, **kwargs):
        return pd.DataFrame()

    def get_disruption_time(self, **kwargs):
        return None

    def get_disruption_shotlist(self, **kwargs):
        return []

    def get_disruption_warning_shotlist(self, **kwargs):
        return []
