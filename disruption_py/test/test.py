import pandas as pd 
import pymysql.cursors
from MySQLdb import _mysql

# import pyodbc
import jaydebeapi
conn = jaydebeapi.connect("com.microsoft.sqlserver.jdbc.SQLServerDriver",
   "jdbc:sqlserver://alcdb2.psfc.mit.edu:1433",
   ["hmturner", "pfcworld"],
   "/home/hmturner/Documents/disruption-warning-db-workflow/cmod/sqljdbc4.jar",)
curs = conn.cursor()
curs.execute("select shot,t_disrupt from disruptions order by shot")
print(curs.fetchall())
curs.close()
# server = 'tcp:alcdb2.psfc.mit.edu,1433'
# database = 'logbook' 
# username = 'hmturner' 
# password = 'pfcworld' 
# cnxn = pyodbc.connect('DRIVER={ODBC Driver 17 for SQL Server};SERVER='+server+';DATABASE='+database+';UID='+username+';PWD='+ password)
# Connect to the database
# connection = pymysql.connect(host='jdbc:sqlserver://alcdb2.psfc.mit.edu:1433',
#                          user='hmturner',
#                              password='pfcworld',
#                              database='logbook',
#                              cursorclass=pymysql.cursors.DictCursor)
# db=_mysql.connect(host="alcdb2.psfc.mit.edu:1433",port=1433,user='hmturner',passwd="pfcworld",db="thangs")

# with connection:
#     with connection.cursor() as cursor:
#         # Create a new record
#         sql = "INSERT INTO `users` (`email`, `password`) VALUES (%s, %s)"
#         cursor.execute(sql, ('webmaster@python.org', 'very-secret'))

#     # connection is not autocommit by default. So you must commit to save
#     # your changes.
#     connection.commit()

#     with connection.cursor() as cursor:
#         # Read a single record
#         sql = "SELECT `id`, `password` FROM `users` WHERE `email`=%s"
#         cursor.execute(sql, ('webmaster@python.org',))
#         result = cursor.fetchone()
#         print(result)

#         db_handle = database(db_name, db_username, db_password, ...
#       'com.microsoft.sqlserver.jdbc.SQLServerDriver', ...
#       ['jdbc:sqlserver://' db_server '.psfc.mit.edu:1433;' ...
#       'database=' db_name]);