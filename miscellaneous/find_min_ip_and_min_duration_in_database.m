db_handle = set_database('logbook');
result = fetch(db_handle, ['select shot, time, ip from disruption_warning ' ...
  'order by shot, time']);
shot = int32(cell2mat(result(:,1)));
time = cell2mat(result(:,2));
ip = cell2mat(result(:,3));
close(db_handle);
clearvars result db_handle;

min_ip = NaN;
min_duration = NaN;

shotlist = unique(shot);
for ishot = 1:length(shotlist);
  indices = find(shot == shotlist(ishot));
  min_ip = min([min_ip max(abs(ip(indices)))]);
  min_duration = min([min_duration max(time(indices))]);
end;

fprintf(1, 'min ip = %6.2f kA\n', min_ip/1.e3);
fprintf(1, 'min duration = %5.1f ms\n', min_duration*1.e3);

clearvars shot time ip shotlist indices;
