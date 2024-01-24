import os
import logging
from typing import List

import pandas as pd
from sqlalchemy import create_engine
# import jaydebeapi
import pyodbc
import threading
        
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
        # self.data_columns = [q[3] for q in self.query("""SELECT *
        # FROM INFORMATION_SCHEMA.COLUMNS
        # WHERE TABLE_NAME = N'disruption_warning'""", use_pandas=True)]

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
            if self._check_protected(query.lower()):
                return 0
        if use_pandas:
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
    
    # Removed until decisions made on updating disruption warnings table
    # def add_shots(self, shots, update=False):
    #     """add_shots to sql table

    #     Parameters
    #     ----------
    #     shots : _type_
    #         _description_
    #     update : bool, optional
    #         _description_, by default False
    #     """
    #     for shot in shots:
    #         self.add_shot(shot=shot, update=update)
            
    # def add_shot(self, shot, update=False):
    #     """
    #     Upload shot to SQL database. Can include shot object if available to avoid redundant computation.
    #     Returns an error if there is at least one row already containing the shot id.
    #     """
    #     if self._check_protected(shot.data.columns):
    #         return False
    #     curr_df = pd.read_sql_query(
    #         f"select * from disruption_warning where shot in {shot.shot_id} order by time", self.engine)
    #     with self.conn.cursor() as curs:
    #         if len(curr_df) == 0:
    #             curs.executemany(
    #                 """insert into disruption_warning values""", shot.data)
    #         elif update:
    #             self._update_shot(shot.data, curr_df, curs)
    #         else:
    #             self.logger.warning("Not updating shot")
    #     return True

    # def _update_shot(shot_data, curr_df, curs):
    #     new_rows = []
    #     for index, row in shot_data.iterrows():
    #         update_row = curr_df[abs(curr_df['time']-row['time']) < TIME_CONST]
    #         if len(update_row) > 1:
    #             # TODO: Change to update all with a new value
    #             raise Exception("Too many matches")
    #         elif len(update_row) == 1:
    #             curs.execute(
    #                 f"""insert into disruption_warning values where dbkey""", update_row, update_row['dbkey'])
    #         else:
    #             new_rows.append(index)
    #     curs.executemany(
    #         f"""insert into disruption_warning values""", curr_df.iloc[new_rows, :])

    # def _check_protected(self, cols_edited):
    #     edits_protected_column = False
    #     for col in self.protected_columns:
    #         if col in cols_edited:
    #             self.logger.error(f"Attempted to edit PROTECTED COLUMN: {col}")
    #             edits_protected_column = True
    #     return not edits_protected_column
        
    
    # # TODO: Protect against injection attacks
    # def remove_shot(self, shot_id):
    #     """ Remove shot from SQL database."""
    #     raise NotImplementedError("Blocking until safety checks are in place")
    #     data_df = pd.read_sql_query(
    #         f'''select * from disruption_warning where shot = {shot_id} order by time''', self.engine)
    #     if len(data_df) == 0:
    #         self.logger.info(f"Shot {shot_id} does not exist in database")
    #         return
    #     with self.conn.cursor() as curs:
    #         curs.execute(
    #             f"delete from disruption_warning where shot = {shot_id}")
    #     return
    
    # def add_column(self, col_name, var_type="TEXT", table="disruption_warning"):
    #     if col_name not in self.data_columns:
    #         self.query(
    #             f"alter table {table} add {col_name} {var_type};", use_pandas=False)
    #         self.data_columns.append(col_name)
    #         return 1
    #     self.logger.info("Column already in database table")
    #     return 0

    # def remove_column(self, col_name, table="disruption_warning"):
    #     if col_name in self.data_columns:
    #         self.query(
    #             f"alter table {table} drop column {col_name};", use_pandas=False)
    #         self.data_columns.remove(col_name)
    #         return 1
    #     self.logger.info("Column not in database table")
    #     return 0

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