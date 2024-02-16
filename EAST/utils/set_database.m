function db_handle = set_database(relational_database);

% The parameter "relational_database" determines which sybase file to open.

if (~exist('relational_database','var'));  % If not specified, then default
  relational_database = 'east_disruption'; % to "east_disruption.sybase_login"
end;

% Get database name, user name, and password from user's sybase file

user = getenv('USER');

fileid = fopen(['/home/ASIPP/' user '/' relational_database '.sybase_login']);
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

% Okay, we should now have everything we need.  Connect to the database.
% (We use the default port for MySql server, which is 3306)

db_handle = database(db_name, db_username, db_password, ...
      'com.mysql.jdbc.Driver', ...
      ['jdbc:mysql://' db_server '.ipp.ac.cn:3306/']);

if (length(db_handle.Message) > 0);
% fprintf(1,['Warning: error accessing %s.ipp.ac.cn\n' ...
%            'Trying 202.127.205.9 instead.\n'], db_server);
  db_handle = database(db_name, db_username, db_password, ...
        'com.mysql.jdbc.Driver', ...
        ['jdbc:mysql://202.127.205.9:3306/']);
end;
