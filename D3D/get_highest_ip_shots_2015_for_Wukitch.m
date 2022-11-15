% Turn off annoying warning messages about java database driver
warning_status = warning;
warning('off')

db = set_database('d3drdb');
result = fetch(db, ['select shot, time, ip from disruption_warning ' ...
  'order by shot, time']);
shot = cell2mat(result(:,1));
time = cell2mat(result(:,2));
ip   = cell2mat(result(:,3));
close(db);

shotlist = unique(shot);
nshots = length(shotlist);

list_max_ip = NaN(nshots,1);
list_min_ip = NaN(nshots,1);

for ishot = 1:nshots;
  indices = find(shot == shotlist(ishot));
  list_max_ip(ishot) = max(ip(indices));
  list_min_ip(ishot) = min(ip(indices));
end;

[list_max_ip, sorted_indices] = sort(list_max_ip, 'descend');
shotlist_max_ip_sorted = shotlist(sorted_indices);

[list_min_ip, sorted_indices] = sort(list_min_ip, 'ascend');
shotlist_min_ip_sorted = shotlist(sorted_indices);

for i=1:5;
  fprintf(1, '%10i  %5.2f MA\n', shotlist_max_ip_sorted(i), ...
    list_max_ip(i)/1e6);
end;

for i=1:5;
  fprintf(1, '%10i  %5.2f MA\n', shotlist_min_ip_sorted(i), ...
    list_min_ip(i)/1e6);
end;

% Turn back on annoying warning messages about java database driver
warning(warning_status);


