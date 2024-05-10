import os
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.constants import CMOD_PROTECTED_COLUMNS


class CModDatabase(ShotDatabase):

    def __init__(self, *args, **kwargs):
        kwargs["protected_columns"] = CMOD_PROTECTED_COLUMNS
        super().__init__(*args, **kwargs)

    @staticmethod
    def default(**kwargs):
        profile = os.path.expanduser("~/logbook.sybase_login")
        with open(profile, "r") as fio:
            config = [line.strip() for line in fio.readlines()]
            # TODO: Account for the optional proxy server on the first line
            _, db_server, db_name, db_user, db_pass = config
        kw = dict(
            driver="ODBC Driver 18 for SQL Server",
            host=db_server,
            port=1433,
            db_name=db_name,
            user=db_user,
            passwd=db_pass,
        )
        kw.update(kwargs)
        db = CModDatabase(**kw)
        return db
