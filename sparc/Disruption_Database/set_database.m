function db_handle = set_database(relational_database);

% The parameter "relational_database" determines which sybase file to open.

if (~exist('relational_database','var')); % If not specified, then
  relational_database = 'logbook';        % default to "logbook.sybase_login"
end;

% Get database name, user name, and password from user's sybase file

user = getenv('USER');

fileid = fopen(['/home/' user '/' relational_database '.sybase_login']);
if (fileid ~= -1); % If sybase file exists in user's directory, read it in
  dummy = fgetl(fileid);
  db_server = fgetl(fileid);
  db_name = fgetl(fileid);
  db_username = fgetl(fileid);
  db_password = fgetl(fileid);
  fclose(fileid);
else; % If sybase file does not exist in user's directory, use default file
  fileid = fopen(['/etc/sybase/' relational_database '.sybase_login']);
  if (fileid ~= -1);
    dummy = fgetl(fileid);
    db_server = fgetl(fileid);
    db_name = fgetl(fileid);
    dummy = fgetl(fileid); % Default sybase file does not contain a username,
    db_username = user;    % so use the login name
    db_password = fgetl(fileid);
    fclose(fileid);
  else;
    fprintf(1,'File "%s.sybase_login" not found.\n', relational_database);
    db_handle = [];
    return;
  end;
end;

% The database connection fails if the outdated servers 'alcdb1' or 'red'
% are specified.  Unfortunately, many people still have 'alcdb1' in their
% sybase file.  So here's a patch to fix the problem for now:

if (strcmpi(db_server, 'alcdb1') || strcmpi(db_server, 'red'));
  db_server = 'alcdb2';
end;

% Define the java class path for the SQL driver, if not already defined.

sql_jdbc_driver = '/home/granetz/java/sqljdbc_4.1/enu/sqljdbc4.jar';

dpath = javaclasspath; % Get current list of dynamic java class path(s)

if sum(strcmpi(sql_jdbc_driver, dpath)) == 0; % If sql driver is not in list,
  javaclasspath(sql_jdbc_driver, dpath);      % add it to list.
end;

% Okay, we should now have everything we need.  Connect to the database.
% (We use the default port for MS sqlserver, which is 1433)

db_handle = database(db_name, db_username, db_password, ...
      'com.microsoft.sqlserver.jdbc.SQLServerDriver', ...
      ['jdbc:sqlserver://' db_server '.psfc.mit.edu:1433;' ...
      'database=' db_name]);
