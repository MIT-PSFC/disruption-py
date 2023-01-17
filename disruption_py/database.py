import os
import logging
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources

import pandas as pd
from pandas.api.types import is_numeric_dtype
import pymysql.cursors
import jaydebeapi
import matplotlib.pyplot as plt

from disruption_py.shots import *
from disruption_py.utils import save_open_plots
import disruption_py.data

# Alter queries for these columns will fail
PROTECTED_COLUMNS = ['dbkey', 'shot', 'time', 'time_until_disrupt', 'ip_error', 'dip_dt', 'beta_p', 'beta_n', 'li', 'n_equal_1_normalized', 'z_error', 'v_z', 'z_times_v_z', 'kappa', 'pressure_peaking', 'H98', 'q0', 'qstar', 'q95', 'v_0', 'v_mid', 'v_edge', 'dn_dt', 'p_rad_slow', 'p_oh_slow', 'p_icrf', 'p_lh', 'radiated_fraction', 'power_supply_railed', 'v_loop_efit', 'r_dd', 'lower_gap', 'upper_gap', 'dbetap_dt', 'dli_dt',
                     'ip', 'zcur', 'n_e', 'dipprog_dt', 'v_loop', 'p_rad', 'p_oh', 'ssep', 'dWmhd_dt', 'dprad_dt', 'v_0_uncalibrated', 'Te_width', 'Greenwald_fraction', 'intentional_disruption', 'Te_width_ECE', 'Wmhd', 'n_over_ncrit', 'n_equal_1_mode', 'Mirnov', 'Mirnov_norm_btor', 'Mirnov_norm_bpol', 'Te_peaking', 'ne_peaking', 'Te_peaking_ECE', 'SXR_peaking', 'kappa_area', 'I_efc', 'SXR', 'H_alpha', 'Prad_peaking_CVA', 'commit_hash']
# [s] Time frame for which an insertion into SQL database becomes an update
TIME_CONST = 1e-6

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class DatabaseHandler:
    """
    Handles grabbing data from MySQL server. 
    """

    def __init__(self, driver, driver_file, host, user, passwd, shot_class=Shot):
        self.user = user
        self.passwd = passwd
        self.driver = driver
        self.driver_file = driver_file
        self.host = host
        self.shot_class = shot_class
        self.conn = jaydebeapi.connect(self.driver,
                                       self.host, [self.user, self.passwd],
                                       self.driver_file)
        self.data_columns = [q[3] for q in self.query("""SELECT *
        FROM INFORMATION_SCHEMA.COLUMNS
        WHERE TABLE_NAME = N'disruption_warning'""", use_pandas=False)]

    def query(self, query, use_pandas=True):
        """ Query SQL database. SELECT queries return a pandas dataframe. """
        if "alter" in query.lower():
            for col in PROTECTED_COLUMNS:
                if col in query.lower():
                    logger.error(f"PROTECTED COLUMN: {col}")
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
            logger.debug(e)
            logger.error("Query failed, returning None")
        curs.close()
        return output

    def add_shots(self, shot_ids, shots=None, update=False):
        shot_ids = tuple([int(shot) for shot in shot_ids])
        if shots is None:
            shots = [Shot(shot_id) for shot_id in shot_ids]
        curr_df = pd.read_sql_query(
            f"select * from disruption_warning where shot in {shot_ids} order by time", self.conn)
        grouped_df = curr_df.groupby(by=["shot"])
        shots_in_db = list(grouped_df.keys())
        # TODO: Paralellize in D3D and EAST subclasses
        for shot in shots:
            self.add_shot(shot._shot_id, shot.data, update=update)

    def add_shot(self, shot_id, shot=None, update=False):
        """ 
        Upload shot to SQL database. Can include shot object if available to avoid redundant computation. 
        Returns an error if there is at least one row already containing the shot id.
        """
        shot_id = int(shot_id)
        if shot is None:
            shot = Shot(shot_id)
        curr_df = pd.read_sql_query(
            f"select * from disruption_warning where shot in {shot_id} order by time", self.conn)
        with self.conn.cursor() as curs:
            if len(curr_df) == 0:
                # self.data.to_sql('disruption_warning', con= self.conn, if_exists='append')
                curs.executemany(
                    """insert into disruption_warning values""", shot.data)
                return
            if update:
                self._update_shot(shot.data, curr_df, curs)
            logger.warning("Not updating shot")

    def _update_shot(shot_data, curr_df, curs):
        new_rows = []
        for index, row in shot.data.iterrows():
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

    # TODO: Protect against injection attacks
    def remove_shot(self, shot_id):
        """ Remove shot from SQL database."""
        raise NotImplementedError("Blocking until safety checks are in place")
        data_df = pd.read_sql_query(
            f'''select * from disruption_warning where shot = {shot_id} order by time''', self.conn)
        if len(data_df) == 0:
            print("Shot does not exist in database")
            return
        with self.conn.cursor() as curs:
            curs.execute(
                f"delete from disruption_warning where shot = {shot_id}")
        return

    def get_shot(self, shot_id):
        shot_id = int(shot_id)
        data_df = pd.read_sql_query(
            f"select * from disruption_warning where shot = {shot_id} order by time", self.conn)
        if self.shot_class == CModShot:
            return CModShot('cmod', shot_id, data=data_df)
        elif self.shot_class == D3DShot:
            # TODO: Better exception
            raise Exception("Invalid shot class for this handler")
        return None

    def get_shots(self, shot_ids=None):
        shot_ids = tuple([int(shot_id) for shot_id in shot_ids])
        if shot_ids is None:
            query = f"select * from disruption_warning order by time"
        else:
            query = f"select * from disruption_warning where shot in {shot_ids} order by time"
        shot_df = pd.read_sql_query(query, self.conn)
        return [Shot('cmod', shot_data['shot'].iloc[0], data=shot_data) for _, shot_data in shot_df.groupby(by=["shot"])]

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

    def get_disruption_shotlist(self):
        """ 
        Get pandas dataframe of all disruption shots and times from the disruption table. Used as a cross-reference to determine whether a given shot is disruptive or not(provides t_disrupt). 
        NOTE: The disruption_warning table contains ONLY a subset of shots in this table
        """
        return self.query('select shot,t_disrupt from disruptions order by shot')

    def get_disruption_table_shotlist(self):
        """ 
        Get pandas dataframe of all shots in the disruption table.
        """
        return self.query('select distinct shot from disruption_warning order by shot')


class D3DHandler(DatabaseHandler):
    def __init__(self, driver, driver_file, host, user, passwd, shot_class=D3DShot):
        super().__init__(driver, driver_file, host +
                         "database=d3drdb", user, passwd, shot_class)
        self.tree_conn = jaydebeapi.connect(self.driver,
                                            host +
                                            "database=code_rundb", [
                                                self.user, self.passwd],
                                            self.driver_file)

    def get_efit_tree(self, shot_id):
        with self.tree_conn.cursor() as curs:
            curs.execute(
                f"select tree from plasmas where shot = {shot_id} and runtag = 'DIS' and deleted = 0 order by idx")
            efit_trees = curs.fetchall()
        efit_tree = efit_trees[-1][0]
        return efit_tree

    def get_shot(self, shot_id, efit_tree=None):
        shot_id = int(shot_id)
        data_df = pd.read_sql_query(
            f"select * from disruption_warning where shot = {shot_id} order by time", self.conn)
        if efit_tree is None:
            efit_tree = self.get_efit_tree(shot_id)
        return D3DShot(shot_id, efit_tree, data=data_df)

    # TODO: Make more efficient
    def get_shots(self, shot_ids, efit_tree=None):
        return [self.get_shot(shot_id, efit_tree) for shot_id in shot_ids]

    def validate_shot(self, shot_id, visualize_differences=False):
        """
        Compare shot data currently in disruption database to what is calculated by the shot object.

        Parameters
        ----------
        shot_id : int or str
            Shot number to validate
        visualize_differences : bool, optional
            Whether to plot the differences between the two dataframes, by default False

        Returns
        -------
        bool
            Whether the shot data is correct according to the disruption database
        """
        shot_id = int(shot_id)
        true_shot = self.get_shot(shot_id)
        if true_shot is None:
            logging.debug("Shot not in database")
            return False
        local_shot = D3DShot(shot_id, self.get_efit_tree(shot_id))
        comparison = pd.DataFrame()
        for col in list(local_shot.data.columns):
            if is_numeric_dtype(local_shot.data[col]):
                comparison[col] = true_shot.data[col] - local_shot.data[col]
        comparison['time'] = true_shot.data['time']
        if not comparison.empty:
            logging.debug("Shot data is not correct")
            for col in comparison.columns:
                if (comparison[col] > 1e-10).any():
                    plt.figure()
                    plt.plot(comparison['time'],
                             comparison[col], label='difference')
                    plt.plot(
                        true_shot.data['time'], true_shot.data[col], label='database')
                    plt.plot(
                        local_shot.data['time'], local_shot.data[col], label='local')
                    plt.title(col)
                    plt.legend()
            save_open_plots('shot_error_graphs.pdf')
            if visualize_differences:
                plt.show()
            return False
        return True


def create_cmod_handler():
    USER = os.getenv('USER')
    # TODO: Catch error if file not found and output helpful error message
    with open(f"/home/{USER}/logbook.sybase_login", "r") as profile:
        content = profile.read().splitlines()[1:]
        db_server = content[0]
        db_name = content[1]
        db_username = content[2]
        assert db_username == USER, f"db_username:{db_username};user:{USER}"
        db_password = content[3]
    with importlib_resources.path(disruption_py.data, "sqljdbc4.jar") as p:
        db_driver_path = str(p)  # Absolute path to jar file
    logging.debug(db_driver_path)
    logging.debug(os.getcwd())
    return DatabaseHandler("com.microsoft.sqlserver.jdbc.SQLServerDriver", db_driver_path, f"jdbc:sqlserver://{db_server}.psfc.mit.edu: 1433", db_username, db_password, shot_class=CmodShot)


def create_d3d_handler():
    USER = os.getenv('USER')
    # TODO: Catch error if file not found and output helpful error message
    with open(f"/home/{USER}/D3DRDB.sybase_login", "r") as profile:
        content = profile.read().splitlines()
        db_username = content[0]
        assert db_username == USER, f"db_username:{db_username};user:{USER}"
        db_password = content[1]
    with importlib_resources.path(disruption_py.data, "sqljdbc4.jar") as p:
        db_driver_path = str(p)
    return D3DHandler("com.microsoft.sqlserver.jdbc.SQLServerDriver", db_driver_path, "jdbc:sqlserver://d3drdb.gat.com:8001;", db_username, db_password, shot_class=D3DShot)


def create_east_handler(self):
    pass


if __name__ == '__main__':
    # test_handler = create_cmod_handler()
    # shot = test_handler.get_shot('1150922001')
    test_handler = create_d3d_handler()
    shot = test_handler.get_shot('175552')
    print(shot.efit_tree_name)
