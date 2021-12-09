function field_names = get_columns(db, table_name);

field_names = fetch(db,['select column_name from ' ...
    'information_schema.columns where table_name = ''' table_name '''']);
