db_retrieve = set_database('d3drdb');
result = fetch(db_retrieve, ['select * from test_disruption_warning ' ...
  'order by shot, time']);
variables = get_columns(db_retrieve, 'test_disruption_warning');
close(db_retrieve);

for ivars = 1:length(variables);
  varname = char(variables(ivars));
  eval([varname '= cell2mat(result(:,ivars));']);
end;
shot = int32(shot); % convert shot numbers to integers
clearvars db_retrieve varname ivars result ;
