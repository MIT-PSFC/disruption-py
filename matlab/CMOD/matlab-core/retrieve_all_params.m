% First, turn off annoying warning messages that tend to be generated when
% using the SQL java driver

warning_status = warning;
warning('off');

db_retrieve = set_database('logbook');
variables = get_columns(db_retrieve, 'disruption_warning');
if strcmp(class(variables), 'table');
  variables = table2cell(variables);
end;
max_charlength = size(char(variables),2);

for ivars = 1:length(variables);
  varname = char(variables(ivars));
  charlength = length(varname);
  strng = [varname blanks(max_charlength - charlength)];
  fprintf(1, 'retrieving %s  (%2i/%2i)\n', strng, ivars, length(variables));
  result_retrieve = fetch(db_retrieve, ['select ' varname ...
    ' from disruption_warning order by shot, time']);
  if strcmp(class(result_retrieve), 'table');
    result_retrieve = table2cell(result_retrieve);
  end;
  eval([varname '= cell2mat(result_retrieve);']);
end;

close(db_retrieve);
shot = int32(shot); % convert shot numbers to integers

% Return warning message status to original setting
if (strcmpi(warning_status(1).state, 'on'));
  warning('on');
end;

clearvars db_retrieve varname ivars result_retrieve max_charlength ...
          charlength strng warning_status;
