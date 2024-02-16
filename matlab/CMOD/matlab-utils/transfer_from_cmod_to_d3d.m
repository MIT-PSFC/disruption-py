db = set_database('logbook');

% Read the D3D data from the C-Mod table.  The shot numbers on D3D
% are <= 6 digits, whereas on C-Mod they are >= 9 digits.  Only 7 columns
% have been populated for each of the D3D shots so far.

values = fetch(db, ['select shot, time, time_until_disrupt, ip, ' ...
  'ip_error, dip_dt, dipprog_dt from disruption_warning where ' ...
  'shot < 1000000 order by shot, time']);
nrecords = size(values, 1);

fields = {'shot', 'time', 'time_until_disrupt', 'ip', 'ip_error', ...
          'dip_dt', 'dipprog_dt'};

% Now insert the data into the d3d database, one record at a time.

for irecord = 1:nrecords
  insert(db, 'disruption_warning_d3d', fields, values(irecord,:));
end;

close(db);
clearvars values nrecords fields irecord db;
