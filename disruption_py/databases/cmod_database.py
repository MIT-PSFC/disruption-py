import os
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.constants import CMOD_PROTECTED_COLUMNS


class CModDatabase(ShotDatabase):

    def __init__(self, driver, host, port, db_name, user, passwd, **kwargs):
        super().__init__(
            driver=driver,
            host=host,
            port=port,
            db_name=db_name,
            user=user,
            passwd=passwd,
            protected_columns=CMOD_PROTECTED_COLUMNS,
            **kwargs
        )

    def default(**kwargs):
        profile = os.path.expanduser("~/logbook.sybase_login")
        with open(profile, "r") as fio:
            db_server, db_name, db_user, db_pass = fio.read().split()
        kw = dict(
            driver="{ODBC Driver 18 for SQL Server}",
            host=db_server,
            port=8001,
            db_name=db_name,
            user=db_user,
            passwd=db_pass,
        )
        kw.update(kwargs)
        db = CModDatabase(**kw)
        return db
