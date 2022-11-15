function db_handle = set_database_cmod(relational_database);

% The parameter "relational_database" determines which SQL database to
% connect to.  Commonly used SQL databases at DIII-D are 'd3drdb' and
% 'code_rundb'.  All SQL databases at DIII-D reside on the server
% 'd3drdb.gat.com'.  Each individual user has a database username and
% password which is contained in a protected file called
% 'D3DRDB.sybase_login' or 'd3drdb.sybase_login' in the user's home
% directory.

% Example of usage:    db = set_database('logbook');

if (~exist('relational_database','var')); % If not specified, then
  relational_database = 'logbook';        % default to C-Mod "logbook" database
end;

% Define java class path to the proper JDBC driver

dpath = getenv('SQLJDBC_JAR'); % on IrisD

javaclasspath('-v0'); % Turn off java status messages
javaclasspath(dpath);

% Get user name and password from user's sybase file

sybase_file = ['/home/' getenv('USER') '/logbook.sybase_login'];
fileid = fopen(sybase_file);
if (fileid == -1); % If sybase file does not exist in user's directory
  fprintf(1,'File "logbook.sybase_login" not found.\n');
  db_handle = [];
  return;
end;

dummy = fgetl(fileid);
db_server = fgetl(fileid);
db_name = fgetl(fileid);
db_username = fgetl(fileid);
db_password = fgetl(fileid);
fclose(fileid);

% Okay, we should now have everything we need.  Connect to the database on
% the Alcator C-Mod database server

db_handle = database(db_name, db_username, db_password, ...
      'com.microsoft.sqlserver.jdbc.SQLServerDriver', ...
      ['jdbc:sqlserver://' db_server '.psfc.mit.edu:1433;' ...
      'database=' db_name]);
