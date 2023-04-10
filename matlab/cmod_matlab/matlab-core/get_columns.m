function field_names = get_columns(db, table_name);

% Determine if the database connection parameter, "db", is a database
% structure or a database handle.  Only the handle is actually needed to
% specify the connection to the database.  By executing this test, this
% routine will work whether or not Matlab's Database Toolbox is available.

%{
if (strcmpi(class(db), 'database'));
  db = db.Handle;  % Only the handle info is needed to access the database
end;
%}

field_names = fetch(db,['select column_name from ' ...
    'information_schema.columns where table_name = ''' table_name '''']);
