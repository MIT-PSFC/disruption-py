db = set_database('logbook');
result = fetch(db, ['select shot from summary where ' ...
  'pulse_length > 0.1 and ipmax >= 0.2e6 and ' ...         
  'shot between 1050101000 and 1051231999 order by shot']);
close(db);
result = table2cell(result);
shotlist = int32(cell2mat(result));
clearvars -except shotlist;
