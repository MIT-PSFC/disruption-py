db = set_database('logbook');
result = fetch(db, 'select * from disruption_warning order by shot, time');
variables = get_columns(db, 'disruption_warning');
close(db);

for ivars = 1:length(variables);
  varname = char(variables(ivars));
  eval([varname '= cell2mat(result(:,ivars));']);
end;
shot = int32(shot); % convert shot numbers to integers
clearvars db varname ivars result ;
