import os
import threading
import pyodbc
from disruption_py.databases import ShotDatabase
from disruption_py.utils.constants import D3D_PROTECTED_COLUMNS


class D3DDatabase(ShotDatabase):

    def __init__(self, driver, host, port, db_name, user, passwd, **kwargs):
        super().__init__(
            driver=driver,
            host=host,
            port=port,
            db_name=db_name,
            user=user,
            passwd=passwd,
            protected_columns=D3D_PROTECTED_COLUMNS,
            **kwargs,
        )
        self._tree_thread_connections = {}
        self.tree_connection_string = self._get_connection_string("code_rundb")

    def default(**kwargs):
        profile = os.path.expanduser("~/D3DRDB.sybase_login")
        with open(profile, "r") as fio:
            db_user, db_pass = fio.read().split()
        kw = dict(
            driver="{FreeTDS}",
            host="d3drdb.gat.com",
            port=8001,
            db_name="d3drdb",
            user=db_user,
            passwd=db_pass,
        )
        kw.update(kwargs)
        db = D3DDatabase(**kw)
        return db

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
        if current_thread not in self._tree_thread_connections:
            self.logger.info(f"Connecting to code_rundb database for thread {current_thread}")
            self._tree_thread_connections[current_thread] = pyodbc.connect(self.tree_connection_string)
        return self._tree_thread_connections[current_thread]

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
