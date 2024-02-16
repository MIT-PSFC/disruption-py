shotlist=dlmread(['/fusion/projects/disruption_warning/matlab_programs/' ...
  'shotlist_rea_blessed.txt']);
shotlist = int32(shotlist);

load('Te_HWHM_data.mat');

db = set_database('d3drdb');

Te_HWHM = NaN(0);
for ishot = 1:length(shotlist);
  result = fetch(db,['select Te_HWHM from disruption_warning ' ...
    'where shot = ' int2str(shotlist(ishot)) ' order by time']);
  Te_HWHM_shot = cell2mat(result);
  Te_HWHM = [Te_HWHM; Te_HWHM_shot];
end;

