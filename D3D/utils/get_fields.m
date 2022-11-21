db = set_database('d3drdb'); % My routine to connect to database

meta = db.getMetaData();
rs = meta.getColumns( {}, {}, 'disruptions', {} );
field = {};
type = {};
nullable ={};

while (rs.next());
  field = [field, rs.getObject('COLUMN_NAME')]; % same as rs.getObject(4)
  type = [type, rs.getObject('TYPE_NAME')];
  nullable = [nullable, rs.getObject('NULLABLE')];
end;

ncols = numel(field);
for i = 1:ncols
  fprintf(1, '%20s  %12s  %5i\n', field{i}, type{i}, nullable{i})
end;

db.close();
rs.close();
