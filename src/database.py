import pandas as pd 
import pymysql.cursors
import jaydebeapi

class DatabaseHandler:
    """
    Handles grabbing data from MySQL server. 
    """
    def __init__(self, driver,driver_file, host, user,passwd):
        self.user = user 
        self.passwd = passwd
        self.driver = driver 
        self.driver_file = driver_file
        self.host = host
        self.conn = jaydebeapi.connect(self.driver,
        self.host,[self.user, self.passwd],
        self.driver_file)
    
    def query(self, query):
        """ Query SQL database """ 
        return pd.read_sql_query(query, self.conn)

    def add_shot(self,shot_id, shot=None):
        """ 
        Upload shot to SQL database. Can include shot object if available to avoid redundant computation. 
        Returns an error if there is at least one row already containing the shot id.
        """
        pass
    def update_shot(self,shot_id, shot=None):
        """

        """
        pass

    def remove_shot(self,shot_id):
        """ Remove shot from SQL database."""
        pass 

    def get_shot(self,shot_id):
        data_df = pd.read_sql_query(f'''select * from disruption_warning where shot = {shot_id}''', self.conn)
        return data_df

    def get_disruption_shotlist():
        """ Get list of all disruption shots and times from the disruption database"""
        pass 
    

def create_cmod_handler():
    return DatabaseHandler("com.microsoft.sqlserver.jdbc.SQLServerDriver","sqljdbc4.jar","jdbc:sqlserver://alcdb2.psfc.mit.edu:1433","hmturner","pfcworld")

def create_d3d_handler():
    pass 

def create_east_handler():
    pass

if __name__ == '__main__':
    test_handler = create_cmod_handler()
    shot = test_handler.get_shot('1150922001')
    print(shot)