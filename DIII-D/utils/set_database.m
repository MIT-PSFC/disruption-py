function db_handle = set_database(relational_database);

% The parameter "relational_database" determines which SQL database to
% connect to.  Commonly used SQL databases at DIII-D are 'd3drdb' and
% 'code_rundb'.  All SQL databases at DIII-D reside on the server
% 'd3drdb.gat.com'.  Each individual user has a database username and
% password which is contained in a protected file called
% 'D3DRDB.sybase_login' or 'd3drdb.sybase_login' in the user's home
% directory.

% Example of usage:    db = set_database('d3drdb');

if (~exist('relational_database','var')); % If not specified, then
  relational_database = 'd3drdb';         % default to "d3drdb" database
end;

% Define java class path to the proper JDBC driver

dpath = getenv('SQLJDBC_JAR'); % on IrisD

javaclasspath('-v0'); % Turn off java status messages
javaclasspath(dpath);

% Get user name and password from user's sybase file

sybase_file = ['/home/' getenv('USER') '/D3DRDB.sybase_login'];
fileid = fopen(sybase_file);
if (fileid == -1); % If sybase file does not exist in user's directory, try
                   % using lowercase file name
  sybase_file = ['/home/' getenv('USER') '/d3drdb.sybase_login'];
  fileid = fopen(sybase_file);
  if (fileid == -1) % Still can't find sybase file
    fprintf(1,'File "D3DRDB.sybase_login" not found.\n');
    db_handle = [];
    return;
  end;
end;

db_username = fgetl(fileid);
db_password = fgetl(fileid);
fclose(fileid);

% Okay, we should now have everything we need.  Connect to the database.
% (We use port 8001 for MS sqlserver.)

db_handle = database(relational_database, db_username, db_password, ...
      'com.microsoft.sqlserver.jdbc.SQLServerDriver', ...
      ['jdbc:sqlserver://d3drdb.gat.com:8001;' ...
      'database=' relational_database]);
