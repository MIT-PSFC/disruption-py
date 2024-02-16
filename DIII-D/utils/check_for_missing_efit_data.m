db1 = set_database('test_db');
result = fetch(db1, ['select shot, time, li ' ...
  'from test_disruption_warning order by shot, time']);
shot = cell2mat(result(:,1));
time = cell2mat(result(:,2));
li = cell2mat(result(:,3));
close(db1);
clearvars result db1;

shotlist = unique(shot);

list_of_missing_efit_data = [];

for ishot = 1:length(shotlist);
  indices = find(shot == shotlist(ishot));
  n_time_slices = length(indices);
  n_efit_data_missing = length(find(isnan(li(indices))));
  frac_efit_data_missing = n_efit_data_missing / n_time_slices;
  if (frac_efit_data_missing > 0.15);
    fprintf(1, '  Missing EFIT data on shot %6i\n', shotlist(ishot));
    list_of_missing_efit_data = [list_of_missing_efit_data, shotlist(ishot)];
  end;
end;

save('list_of_missing_efit_data.mat', 'list_of_missing_efit_data');
fileid = fopen('list_of_missing_efit_data.txt', 'w');
  fprintf(fileid, '%6i\n', list_of_missing_efit_data);
fclose(fileid);
