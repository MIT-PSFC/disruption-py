db = set_database('logbook');
result = fetch(db,['select shot, time, time_until_disrupt ' ...
  'from disruption_warning_d3d order by shot, time']);
close(db);
shots=cell2mat(result(:,1));
times=cell2mat(result(:,2));
tdis =cell2mat(result(:,3));
clearvars result;

shotlist = unique(shots);
end_of_shot = NaN(size(shotlist));
disrupt_time = NaN(size(shotlist));

for ishot = 1:length(shotlist);
  indices = find(shots == shotlist(ishot));
  [end_of_shot(ishot), indx] = max(times(indices));
  if (~isnan(tdis(indices(indx))));
    disrupt_time(ishot) = end_of_shot(ishot);
  end;
end;

save('d3d_shot_duration_disrupt.mat', 'shotlist', 'end_of_shot', ...
  'disrupt_time');
