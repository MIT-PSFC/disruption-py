import os
import logging
from typing import List

import pandas as pd
import numpy as np
from sqlalchemy import create_engine
# import jaydebeapi
import pyodbc
import threading

from disruption_py.utils.constants import TIME_CONST
        
class ShotDatabase:
    """
    Handles grabbing data from MySQL server.
    """
    logger = logging.getLogger('disruption_py')        
        
    def __init__(self, driver, host, db_name, user, passwd, protected_columns=[], **kwargs):
        self.user = user
        self.passwd = passwd
        self.driver = driver
        self.host = host
        self.db_name = db_name
        self.protected_columns = protected_columns
        self.connection_string = (
            f"DRIVER={self.driver};"
            f"SERVER={self.host};"
            f"DATABASE={self.db_name};"
            f"UID={self.user};"
            f"PWD={self.passwd};"
            "TrustServerCertificate=yes;"
            "Connection Timeout=60"
        )
        self._thread_connections = {}
        self.logger.info("Database initialized")
        self.engine = create_engine(f"mssql+pyodbc:///?odbc_connect={self.connection_string}")

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
            self._thread_connections[current_thread] = pyodbc.connect(self.connection_string)
        return self._thread_connections[current_thread]
    
    def query(self, query : str, use_pandas=True):
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
            print(e)
            self.logger.debug(e)
            self.logger.error("Query failed, returning None")
        curs.close()
        return output
            
    def add_shot_data(self, shot_id : int, shot_data : pd.DataFrame, update=False, override_update_columns=None):
        """
        Upload shot to SQL database. Can include shot object if available to avoid redundant computation.
        Returns an error if there is at least one row already containing the shot id.
        """
        
        curr_df = pd.read_sql_query(
            f"select * from disruption_warning where shot in {shot_id} order by time", self.engine)
        
        if len(curr_df == 0):
            return self._insert_shot_data(curr_df, shot_data)
        elif len(curr_df) == len(shot_data) and ((curr_df['time']  - shot_data['time']).abs() < TIME_CONST).all():
            return self._update_shot_data(curr_df, shot_data, update, override_update_columns)
        
        self.logger.error("Invalid timebase for data output")
        return False

    def _insert_shot_data(self, curr_df : pd.DataFrame, shot_data : pd.DataFrame):
        matching_columns_shot_data = pd.DataFrame()
        for column_name in curr_df.columns:
            if column_name in shot_data.columns:
                matching_columns_shot_data[column_name] = shot_data[column_name]
            else:
                matching_columns_shot_data[column_name] = np.nan
        with self.conn.cursor() as curs:
            data_tuples = list(matching_columns_shot_data.itertuples(index=False, name=None))
            curs.executemany(
                """insert into disruption_warning values""", data_tuples)
        return True
    
    def _update_shot_data(self, curr_df : pd.DataFrame, shot_data : pd.DataFrame, update : bool, override_update_columns=None):
        NON_UPDATABLE_COLUMNS = ["shot", "time"]
        override_update_columns = override_update_columns or []
        
        update_columns_shot_data = pd.DataFrame()
        for column_name in curr_df.columns:
            if (
                column_name in NON_UPDATABLE_COLUMNS or 
                (column_name in self.protected_columns and column_name not in override_update_columns)
            ):
                continue
            
            if (
                column_name in shot_data.columns and
                not shot_data[column_name].isna().all() and
                (update or curr_df[column_name].isna().all())
            ):
                update_columns_shot_data[column_name] = shot_data[column_name]
        with self.conn.cursor() as curs:
            for index, row in update_columns_shot_data.iterrows():     
                update_columns = update_columns_shot_data.columns.difference(NON_UPDATABLE_COLUMNS)
                sql_command = f"UPDATE disruption_warning SET {', '.join([f'{col} = ?' for col in update_columns])} WHERE time = ?;"
                curs.execute(sql_command, row + curr_df['time'][index])
        return True
        
    
    def remove_shot_data(self, shot_id):
        """Remove shot from SQL database."""
        data_df = pd.read_sql_query(
            f'''select * from disruption_warning where shot = {shot_id} order by time''', self.engine)
        if len(data_df) == 0:
            self.logger.info(f"Shot {shot_id} does not exist in database")
            return False
        with self.conn.cursor() as curs:
            curs.execute(
                f"delete from disruption_warning where shot = {shot_id}")
        return True
    
    def add_column(self, col_name, var_type="TEXT", table="disruption_warning"):
        try:
            self.query(f"alter table {table} add {col_name} {var_type};", use_pandas=False)
            return True
        except Exception as e:
            self.logger.error(f"Failed to add column with error {e}")
            return False

    def remove_column(self, col_name, table="disruption_warning"):
        try:
            self.query(f"alter table {table} drop column {col_name};", use_pandas=False)
            return True
        except Exception as e:
            self.logger.error(f"Failed to drop column with error {e}")
            return False

    def get_shot_data(self, shot_ids : List[int], cols : List[str]=["*"], sql_table="disruption_warning"):
        """get_shot_data retrieves columns from sql data for given shot_ids

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
        shot_ids = ','.join([str(shot_id) for shot_id in shot_ids])
        selected_cols = f"{cols[0]}"
        if len(cols) > 1:
            selected_cols += ''.join([f", {col}" for col in cols[1:]])
        if shot_ids is None:
            query = f"select {selected_cols} from {sql_table} order by time"
        else:
            query = f"select {selected_cols} from {sql_table} where shot in ({shot_ids}) order by time"
        shot_df = pd.read_sql_query(query, self.engine)
        return shot_df

    def get_disruption_time(self, shot_id):
        """
        Get disruption time for shot_id or None if there was no disruption.
        """
        with self.conn.cursor() as curs:
            curs.execute(
                f"select t_disrupt from disruptions where shot = {shot_id}")
            t_disrupt = curs.fetchall()
        if len(t_disrupt) == 0:
            return None
        t_disrupt = t_disrupt[0][0]
        return t_disrupt

    def get_disruption_shotlist(self):
        """ 
        Get pandas dataframe of all disruptive shots and times from the disruption table. Can be sed as a cross-reference to determine whether a given shot is disruptive or not (all shots in this table are disruptive) and contain a t_disrupt.  
        """
        return self.query('select distinct shot from disruptions order by shot')

    def get_disruption_warning_shotlist(self):
        """ 
        Get pandas dataframe of all shots in the disruption_warning table. NOTE: The disruption_warning table contains ONLY a subset of shots in this table
        """
        return self.query('select distinct shot from disruption_warning order by shot')