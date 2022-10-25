import pandas as pd 
import pymysql.cursors
import jaydebeapi
from plasma import *

# Alter queries for these columns will fail 
PROTECTED_COLUMNS = ['dbkey', 'shot', 'time', 'time_until_disrupt', 'ip_error', 'dip_dt', 'beta_p', 'beta_n', 'li', 'n_equal_1_normalized', 'z_error', 'v_z', 'z_times_v_z', 'kappa', 'pressure_peaking', 'H98', 'q0', 'qstar', 'q95', 'v_0', 'v_mid', 'v_edge', 'dn_dt', 'p_rad_slow', 'p_oh_slow', 'p_icrf', 'p_lh', 'radiated_fraction', 'power_supply_railed', 'v_loop_efit', 'r_dd', 'lower_gap', 'upper_gap', 'dbetap_dt', 'dli_dt', 'ip', 'zcur', 'n_e', 'dipprog_dt', 'v_loop', 'p_rad', 'p_oh', 'ssep', 'dWmhd_dt', 'dprad_dt', 'v_0_uncalibrated', 'Te_width', 'Greenwald_fraction', 'intentional_disruption', 'Te_width_ECE', 'Wmhd', 'n_over_ncrit', 'n_equal_1_mode', 'Mirnov', 'Mirnov_norm_btor', 'Mirnov_norm_bpol', 'Te_peaking', 'ne_peaking', 'Te_peaking_ECE', 'SXR_peaking', 'kappa_area', 'I_efc', 'SXR', 'H_alpha', 'Prad_peaking_CVA', 'commit_hash']

class DatabaseHandler:
    """
    Handles grabbing data from MySQL server. 
    """
    def __init__(self, driver,driver_file, host, user,passwd,shot_class=Shot):
        self.user = user 
        self.passwd = passwd
        self.driver = driver 
        self.driver_file = driver_file
        self.host = host
        self.shot_class = Shot
        self.conn = jaydebeapi.connect(self.driver,
        self.host,[self.user, self.passwd],
        self.driver_file)
        self.data_columns = [ q[3] for q in self.query("""SELECT *
        FROM INFORMATION_SCHEMA.COLUMNS
        WHERE TABLE_NAME = N'disruption_warning'""", use_pandas=False)]
    
    def query(self, query,use_pandas=True):
        """ Query SQL database. SELECT queries return a pandas dataframe. """
        if "alter" in query.lower():
            for col in PROTECTED_COLUMNS:
                if col in query.lower():
                    print("ERROR: PROTECTED COLUMN")
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
            print("ERROR: Query failed, returning None")
        curs.close()
        return output

    def add_shot(self,shot_id, shot=None):
        """ 
        Upload shot to SQL database. Can include shot object if available to avoid redundant computation. 
        Returns an error if there is at least one row already containing the shot id.
        """
        pass

    def update_shot(self,shot_id, shot=None):
        """ """
        pass

    def remove_shot(self,shot_id):
        """ Remove shot from SQL database."""
        pass 

    def get_shot(self,shot_id):
        shot_id = int(shot_id)
        data_df = pd.read_sql_query(f'''select * from disruption_warning where shot = {shot_id}''', self.conn)
        if self.shot_class == Shot:
            return Shot('cmod',shot_id,data=data_df)
        return None

    def add_column(self, col_name, var_type="TEXT", table="disruption_warning"):
        if col_name not in self.data_columns:
            self.query(f"alter table {table} add {col_name} {var_type};",use_pandas=False)
            self.data_columns.append(col_name)
            return 1 
        print("Column already in database table")
        return 0
    
    def remove_column(self, col_name, table="disruption_warning"):
        if col_name in self.data_columns:
            self.query(f"alter table {table} drop column {col_name};",use_pandas=False)
            self.data_columns.remove(col_name)
            return 1 
        print("Column not in database table")
        return 0

    def get_disruption_shotlist(self):
        """ 
        Get pandas dataframe of all disruption shots and times from the disruption table. Used as a cross-reference to determine whether a given shot is disruptive or not(provides t_disrupt). 
        NOTE: The disruption_warning table contains ONLY a subset of shots in this table
        """
        return self.query('select shot,t_disrupt from disruptions order by shot')
    
    @classmethod
    def create_cmod_handler(self):
        return DatabaseHandler("com.microsoft.sqlserver.jdbc.SQLServerDriver","../src/sqljdbc4.jar","jdbc:sqlserver://alcdb2.psfc.mit.edu:1433","hmturner","pfcworld")

    @classmethod
    def create_d3d_handler(self):
        pass 

    @classmethod
    def create_east_handler(self):
        pass

if __name__ == '__main__':
    test_handler = create_cmod_handler()
    shot = test_handler.get_shot('1150922001')
    print(shot)