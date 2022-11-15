warning('off');

db1 = set_database('d3drdb');
db2 = set_database('test_db');

fields1 = get_columns(db1, 'disruption_warning');
fields2 = get_columns(db2, 'test_disruption_warning');

if (length(fields1) ~= length(fields2));
  fprintf(1, 'Tables have different number of fields\n');
  close(db1);
  close(db2);
  return;
end;

for i=1:length(fields1);
  if (~strcmp(fields1(i), fields2(i)));
    fprintf(1,'%i  variable lists do not match\n', i);
    close(db1);
    close(db2);
    return;
  end;
end;

no_key_indices = find(~strcmp(fields1, {'dbkey'}));

record_keys = cell2mat(fetch(db2, ['select dbkey from ' ...
  'test_disruption_warning order by shot, time']));
nrecords = length(record_keys);

for i=1:nrecords;
  fprintf(1, 'Transferring record # %6i out of %6i (%6.2f)\n', ...
    i, nrecords, i/nrecords*100);
  values = fetch(db2,['select * from test_disruption_warning where ' ...
    'dbkey = ' int2str(record_keys(i)) ]);
  insert(db1, 'disruption_warning', fields2(no_key_indices), ...
    values(no_key_indices));
end;

close(db1);
close(db2);
warning('on');
