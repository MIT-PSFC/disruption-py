import os
import logging

import pandas as pd
import jaydebeapi

logger = logging.getLogger('disruption_py')
TIME_CONST = 1e-6

class ShotDatabase:
    """
    Handles grabbing data from MySQL server.
    """

    def __init__(self, driver, driver_file, host, user, passwd, protected_columns=[], **kwargs):
        self.user = user
        self.passwd = passwd
        self.driver = driver
        self.driver_file = driver_file
        self.host = host
        self.protected_columns = protected_columns
        self.conn = jaydebeapi.connect(self.driver,
                                       self.host, [self.user, self.passwd],
                                       self.driver_file)
        self.data_columns = [q[3] for q in self.query("""SELECT *
        FROM INFORMATION_SCHEMA.COLUMNS
        WHERE TABLE_NAME = N'disruption_warning'""", use_pandas=False)]

    def query(self, query, use_pandas=True):
        """ Query SQL database. SELECT queries return a pandas dataframe. """
        if "alter" in query.lower():
            if self._check_protected(query.lower()):
                return 0
        if use_pandas:
            return pd.read_sql_query(query, self.conn)
        curs = self.conn.cursor()
        output = None
        try:
            curs.execute(query)
            if "select" in query.lower():
                output = curs.fetchall()
        except jaydebeapi.DatabaseError as e:
            print(e)
            logger.debug(e)
            logger.error("Query failed, returning None")
        curs.close()
        return output
    
    def add_shots(self, shots, update=False):
        for shot in shots:
            self.add_shot(shot=shot, update=update)
            
    def add_shot(self, shot, update=False):
        """
        Upload shot to SQL database. Can include shot object if available to avoid redundant computation.
        Returns an error if there is at least one row already containing the shot id.
        """
        if self._check_protected(shot.data.columns):
            return False
        curr_df = pd.read_sql_query(
            f"select * from disruption_warning where shot in {shot.shot_id} order by time", self.conn)
        with self.conn.cursor() as curs:
            if len(curr_df) == 0:
                curs.executemany(
                    """insert into disruption_warning values""", shot.data)
            elif update:
                self._update_shot(shot.data, curr_df, curs)
            else:
                logger.warning("Not updating shot")
        return True

    def _update_shot(shot_data, curr_df, curs):
        new_rows = []
        for index, row in shot_data.iterrows():
            update_row = curr_df[abs(curr_df['time']-row['time']) < TIME_CONST]
            if len(update_row) > 1:
                # TODO: Change to update all with a new value
                raise Exception("Too many matches")
            elif len(update_row) == 1:
                curs.execute(
                    f"""insert into disruption_warning values where dbkey""", update_row, update_row['dbkey'])
            else:
                new_rows.append(index)
        curs.executemany(
            f"""insert into disruption_warning values""", curr_df.iloc[new_rows, :])

    def _check_protected(self, cols_edited):
        edits_protected_column = False
        for col in self.protected_columns:
            if col in cols_edited:
                logger.error(f"Attempted to edit PROTECTED COLUMN: {col}")
                edits_protected_column = True
        return not edits_protected_column
        
    
    # TODO: Protect against injection attacks
    def remove_shot(self, shot_id):
        """ Remove shot from SQL database."""
        raise NotImplementedError("Blocking until safety checks are in place")
        data_df = pd.read_sql_query(
            f'''select * from disruption_warning where shot = {shot_id} order by time''', self.conn)
        if len(data_df) == 0:
            logger.info(f"Shot {shot_id} does not exist in database")
            return
        with self.conn.cursor() as curs:
            curs.execute(
                f"delete from disruption_warning where shot = {shot_id}")
        return

    def get_shot_data(self, shot_ids=None, cols=["*"], sql_table="disruption_warning"):
        shot_ids = ','.join([str(shot_id) for shot_id in shot_ids])
        selected_cols = f"{cols[0]}"
        if len(cols) > 1:
            selected_cols += ''.join([f", {col}" for col in cols[1:]])
        if shot_ids is None:
            query = f"select {selected_cols} from {sql_table} order by time"
        else:
            query = f"select {selected_cols} from {sql_table} where shot in ({shot_ids}) order by time"
        shot_df = pd.read_sql_query(query, self.conn)
        return shot_df

    def add_column(self, col_name, var_type="TEXT", table="disruption_warning"):
        if col_name not in self.data_columns:
            self.query(
                f"alter table {table} add {col_name} {var_type};", use_pandas=False)
            self.data_columns.append(col_name)
            return 1
        logger.info("Column already in database table")
        return 0

    def remove_column(self, col_name, table="disruption_warning"):
        if col_name in self.data_columns:
            self.query(
                f"alter table {table} drop column {col_name};", use_pandas=False)
            self.data_columns.remove(col_name)
            return 1
        logger.info("Column not in database table")
        return 0

    def get_disruption_time(self, shot_id):
        with self.conn.cursor() as curs:
            curs.execute(
                f"select t_disrupt from disruptions where shot = {shot_id}")
            t_disrupt = curs.fetchall()
        if len(t_disrupt) == 0:
            return None
        t_disrupt = t_disrupt[0][0]
        return t_disrupt

    def get_disruption_time(self, shot_id):
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