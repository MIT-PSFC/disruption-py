import os
import logging
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources

import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np
import pymysql.cursors
import jaydebeapi
import matplotlib.pyplot as plt

from disruption_py.shots import *
from disruption_py.utils import save_open_plots
import disruption_py.data

# Alter queries for these columns will fail
CMOD_PROTECTED_COLUMNS = ['dbkey', 'shot', 'time', 'time_until_disrupt', 'ip_error', 'dip_dt', 'beta_p', 'beta_n', 'li', 'n_equal_1_normalized', 'z_error', 'v_z', 'z_times_v_z', 'kappa', 'pressure_peaking', 'H98', 'q0', 'qstar', 'q95', 'dn_dt', 'p_rad_slow', 'p_oh_slow', 'p_icrf', 'p_lh', 'radiated_fraction', 'power_supply_railed', 'v_loop_efit', 'lower_gap', 'upper_gap', 'dbetap_dt', 'dli_dt',
                          'ip', 'zcur', 'n_e', 'dipprog_dt', 'v_loop', 'p_rad', 'p_oh', 'ssep', 'dWmhd_dt', 'dprad_dt', 'Te_width', 'Greenwald_fraction', 'intentional_disruption', 'Te_width_ECE', 'Wmhd', 'n_over_ncrit', 'n_equal_1_mode', 'Mirnov', 'Mirnov_norm_btor', 'Mirnov_norm_bpol', 'Te_peaking', 'ne_peaking', 'Te_peaking_ECE', 'SXR_peaking', 'kappa_area', 'I_efc', 'SXR', 'H_alpha', 'Prad_peaking_CVA', 'commit_hash']
D3D_PROTECTED_COLUMNS = []
# [s] Time frame for which an insertion into SQL database becomes an update
TIME_CONST = 1e-6

logger = logging.getLogger('disruption_py')


class DatabaseHandler:
    """
    Handles grabbing data from MySQL server.
    """

    def __init__(self, driver, driver_file, host, user, passwd, shot_class=Shot, protected_columns=[]):
        self.user = user
        self.passwd = passwd
        self.driver = driver
        self.driver_file = driver_file
        self.host = host
        self.shot_class = shot_class
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
            for col in self.protected_columns:
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
            print(e)
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
            logger.info(f"Shot {shot_id} does not exist in database")
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
            logger.info(f"Grabbed shot: {shot_id}")
            return CModShot('cmod', shot_id, data=data_df)
        elif self.shot_class == D3DShot:
            # TODO: Better exception
            raise Exception("Invalid shot class for this handler")
        return None

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

    def get_shots(self, shot_ids=None):
        shot_df = self.get_shot_data(shot_ids)
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

class CModHandler(DatabaseHandler):
    def get_shot(self, shot_id, efit_tree=None):
        shot_id = int(shot_id)
        data_df = pd.read_sql_query(
            f"select * from disruption_warning where shot = {shot_id} order by time", self.conn)
        if efit_tree is None:
            efit_tree = 'analysis'
        return CModShot(shot_id, efit_tree_name=efit_tree, data=data_df,disruption_time=self.get_disruption_time(shot_id))

    # TODO: Make more efficient
    def get_shots(self, shot_ids, efit_tree=None):
        return [self.get_shot(shot_id, efit_tree) for shot_id in shot_ids] 

class D3DHandler(DatabaseHandler):
    def __init__(self, driver, driver_file, host, user, passwd, shot_class=D3DShot):
        super().__init__(driver, driver_file, host +
                         "database=d3drdb", user, passwd, shot_class, protected_columns=D3D_PROTECTED_COLUMNS)
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
        if len(efit_trees) == 0:
            efit_trees = [('EFIT01',)]
            # with self.tree_conn.cursor() as curs:
                # curs.execute(f"select tree from plasmas where shot = {shot_id} and deleted = 0 order by idx")
                # efit_trees = curs.fetchall()
        efit_tree = efit_trees[-1][0]
        return efit_tree

    # TODO: Implement
    def get_efit_trees(self, shot_ids):
        pass

    def get_shot(self, shot_id, efit_tree=None):
        shot_id = int(shot_id)
        data_df = pd.read_sql_query(
            f"select * from disruption_warning where shot = {shot_id} order by time", self.conn)
        if efit_tree is None:
            efit_tree = self.get_efit_tree(shot_id)
            if efit_tree is None:
                logging.info(
                    f"Shot {shot_id} has no disruptions group efit run")
                return None
        return D3DShot(shot_id, efit_tree, data=data_df,disruption_time=self.get_disruption_time(shot_id))

    # TODO: Make more efficient
    def get_shots(self, shot_ids, efit_tree=None):
        return [self.get_shot(shot_id, efit_tree) for shot_id in shot_ids]

    def validate_shot(self, shot_id, visualize_differences=False, output_dir='./'):
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
        try:
            true_shot = self.get_shot(shot_id)
            if true_shot is None:
                logging.debug("Shot not in database")
                return False, "Shot not in database"
        except Exception as e:
            logging.warning(f"Failed to load shot {shot_id} from SQL database")
            logging.debug(e)
            return False, "Failed to load shot"
        local_shot = D3DShot(shot_id, self.get_efit_tree(
            shot_id), disruption_time=self.get_disruption_time(shot_id))
        comparison = pd.DataFrame()
        for col in list(local_shot.data.columns):
            if col in list(true_shot.data.columns):
                if is_numeric_dtype(local_shot.data[col]):
                    comparison[col] = np.abs((true_shot.data[col] - local_shot.data[col])/true_shot.data[col])
        validation = False
        if not comparison.empty:
            logging.debug("Shot data is not correct")
            for col in comparison.columns:
                if (comparison[col] > 1e-4).any():
                    fig, axs = plt.subplots(2, 1, sharex=True)
                    axs[0].plot(
                        true_shot.data['time'], true_shot.data[col], label='database')
                    axs[0].plot(
                        local_shot.data['time'], local_shot.data[col], label='local', linestyle="dotted")
                    axs[1].plot(local_shot.data['time'], comparison[col])
                    axs[0].set_title(f"Shot {shot_id}:{col}")
                    axs[1].set_title("Normalized Error")
                    axs[0].legend()
            save_open_plots(output_dir + f"shot_{shot_id}_error_graphs.pdf")
            plt.close('all')
            if visualize_differences:
                plt.show()
        return validation, comparison


def create_cmod_handler():
    USER = os.getenv('USER')
    # TODO: Catch error if file not found and output helpful error message
    with open(f"/home/{USER}/logbook.sybase_login", "r") as profile:
        content = profile.read().splitlines()[1:]
        db_server = content[0]
        db_name = content[1]
        db_username = content[2]
        #assert db_username == USER, f"db_username:{db_username};user:{USER}"
        db_password = content[3]
    with importlib_resources.path(disruption_py.data, "sqljdbc4.jar") as p:
        db_driver_path = str(p)  # Absolute path to jar file
    logging.debug(db_driver_path)
    logging.debug(os.getcwd())
    return CModHandler("com.microsoft.sqlserver.jdbc.SQLServerDriver", db_driver_path, f"jdbc:sqlserver://{db_server}.psfc.mit.edu: 1433", db_username, db_password, shot_class=CModShot)


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
    shot = test_handler.get_shot('180808')
    print(shot.efit_tree_name)
