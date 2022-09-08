db_handle = set_database('logbook');
result = fetch(db_handle, ['select shot, ipmax, pulse_length from ' ...
  'summary where shot between 1160101000 and 1161231999 order by shot']);
shot = int32(cell2mat(result(:,1)));
ipmax = cell2mat(result(:,2));
pulse_length = cell2mat(result(:,3));
close(db_handle);
clearvars result db_handle;

indices = find(abs(ipmax) >= 100e3 & pulse_length > 0.100);
shotlist_2016 = shot(indices);

clearvars shot ipmax pulse_length indices;
