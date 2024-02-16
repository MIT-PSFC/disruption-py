function table_names = get_tables(db);

table_names = fetch(db,['select table_name from ' ...
    'information_schema.tables where table_type = ''base table''']);
