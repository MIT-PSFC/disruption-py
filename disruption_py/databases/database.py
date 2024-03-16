import os
import logging
from typing import List

import pandas as pd
import numpy as np
from sqlalchemy import create_engine
# import jaydebeapi
import pyodbc
import threading

from disruption_py.utils.constants import BASE_PROTECTED_COLUMNS, TIME_CONST
        
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
            
    def add_shot_data(self, shot_id : int, shot_data : pd.DataFrame, update=False, override_columns:List[str]=None, table_name="disruption_warning"):
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
        curr_df = pd.read_sql_query(
            f"select * from {table_name} where shot={shot_id} order by time", self.engine)
        
        if len(curr_df) == 0:
            return self._insert_shot_data(curr_df, shot_data, table_name=table_name)
        elif len(curr_df) == len(shot_data) and ((curr_df['time']  - shot_data['time']).abs() < TIME_CONST).all():
            return self._update_shot_data(shot_id, curr_df, shot_data, update, override_columns, table_name=table_name)
        
        self.logger.error("Invalid timebase for data output")
        return False

    def _insert_shot_data(self, curr_df : pd.DataFrame, shot_data : pd.DataFrame, table_name="disruption_warning"):
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
            data_tuples = list(matching_columns_shot_data.itertuples(index=False, name=None))
            curs.executemany(f"insert into {table_name} ({sql_column_names}) values {parameter_markers}", data_tuples)
        return True
    
    def _update_shot_data(self, shot_id : int, curr_df : pd.DataFrame, shot_data : pd.DataFrame, update : bool, override_columns:List[str]=None, table_name="disruption_warning"):
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
            List of columns that can should still be updated. Update must be true for existing values in the 
            columns to be changed. Default value is [].
        table_name : str
            Name of the table for data insert or update. Default value is "disruption_warning".
        """
        override_columns = override_columns or []
        
        update_columns_shot_data = pd.DataFrame()
        for column_name in curr_df.columns:
            if (
                column_name in BASE_PROTECTED_COLUMNS or 
                (column_name in self.protected_columns and column_name not in override_columns)
            ):
                continue
            
            if (
                column_name in shot_data.columns and
                not shot_data[column_name].isna().all() and
                (update or curr_df[column_name].isna().all())
            ):
                update_columns_shot_data[column_name] = shot_data[column_name]
        with self.conn.cursor() as curs:
            for index, row in enumerate(update_columns_shot_data.itertuples(index=False, name=None)):     
                update_column_names = list(update_columns_shot_data.columns)
                sql_set_string = ', '.join([f'{col} = ?' for col in update_column_names])
                sql_command = f"UPDATE {table_name} SET {sql_set_string} WHERE time = ? AND shot = ?;"
                curs.execute(sql_command, row + (curr_df['time'][index], str(shot_id)))
        return True
    
    def _get_identity_column_names(self, table_name="disruption_warning"):
        """Get which column names are identity columns in table."""
        with self.conn.cursor() as curs:
            query = f"""\
            SELECT c.name AS ColumnName 
            FROM sys.columns c
            INNER JOIN sys.tables t ON c.object_id = t.object_id
            LEFT JOIN sys.identity_columns ic ON ic.object_id = c.object_id AND ic.column_id = c.column_id
            WHERE t.name = '{table_name}' AND ic.object_id IS NOT NULL
            """
            curs.execute(query)
            return [row[0] for row in curs.fetchall()]
    
    def remove_shot_data(self, shot_id, table_name="disruption_warning"):
        """Remove shot from SQL table."""
        if table_name == "disruption_warning":
            raise ValueError("Please do not delete from the disruption_warning database")
        data_df = pd.read_sql_query(
            f'''select * from {table_name} where shot = {shot_id} order by time''', self.engine)
        if len(data_df) == 0:
            self.logger.info(f"Shot {shot_id} does not exist in database")
            return False
        with self.conn.cursor() as curs:
            curs.execute(
                f"delete from {table_name} where shot = {shot_id}")
        return True
    
    def add_column(self, col_name, var_type="TEXT", table_name="disruption_warning"):
        """Add column to SQL table without filling in data for column."""
        try:
            self.query(f"alter table {table_name} add {col_name} {var_type};", use_pandas=False)
            return True
        except Exception as e:
            self.logger.error(f"Failed to add column {col_name} with error {e}")
            return False

    def remove_column(self, col_name, table_name="disruption_warning"):
        """Remove column from SQL table"""
        if col_name in self.protected_columns:
            self.logger.error(f"Failed to drop protected column {col_name}")
            return False
        try:
            self.query(f"alter table {table_name} drop column {col_name};", use_pandas=False)
            return True
        except Exception as e:
            self.logger.error(f"Failed to drop column {col_name} with error {e}")
            return False

    def get_shots_data(self, shot_ids : List[int], cols : List[str]=["*"], sql_table="disruption_warning"):
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