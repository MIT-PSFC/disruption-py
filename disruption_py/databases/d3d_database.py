import os
import logging
import threading
import pyodbc
from disruption_py.databases import ShotDatabase
from disruption_py.utils.constants import D3D_PROTECTED_COLUMNS

class D3DDatabase(ShotDatabase):
    logger = logging.getLogger('disruption_py')

    def __init__(self, driver, host, db_name, user, passwd, **kwargs):
        super().__init__(driver, host, db_name, user, passwd, protected_columns=D3D_PROTECTED_COLUMNS, **kwargs)
        self.tree_connection_string = self._get_connection_string("code_rundb")
  
    def default(**kwargs):
        USER = os.getenv('USER')
        # TODO: Catch error if file not found and output helpful error message
        with open(f"/home/{USER}/D3DRDB.sybase_login", "r") as profile:
            content = profile.read().splitlines()
            db_username = content[0]
            assert db_username == USER, f"db_username:{db_username};user:{USER}"
            db_password = content[1]
        return D3DDatabase(
            driver="FreeTDS", 
            host= "d3drdb.gat.com:8001",
            db_name= "D3DRDB",
            user=db_username, 
            passwd=db_password,
            **kwargs
        )
  
    @property
    def tree_conn(self):
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
