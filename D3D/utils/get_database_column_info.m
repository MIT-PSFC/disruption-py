db = set_database('d3drdb');
result=fetch(db,['select * from information_schema.columns ' ...
  'where table_name=''disruption_warning'' ']);

result(:,[4,6,7,8]) % print table information

column_name = result(:,4);
column_default = result(:,6);
is_nullable = result(:,7);
data_type = result(:,8);
