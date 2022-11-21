db = set_database('test_db');
shotlist = cell2mat(fetch(db, ['select distinct(shot) from ' ...
  'test_disruption_warning order by shot']));
nshots = length(shotlist);

maxcount = 0;
for ishot = 1:nshots;
  shot = shotlist(ishot);
  fprintf(1, 'Processing shot %6i [%4i/%4i  (%6.2f%%)]\n', shot, ishot, ...
    nshots, ishot/nshots*100);
  time = cell2mat(fetch(db, ['select time from test_disruption_warning ' ...
    'where shot = ' num2str(shot) ' order by time']));
  count = length(find(diff(time) <= 10e-6));
  maxcount = max([maxcount count]);
  if (count > 0);
    fprintf(1, '  Shot %6i had %4i duplicate times\n', shot, count);
  end;
end;

close(db);
fprintf(1, 'maxcount = %i\n', maxcount);
