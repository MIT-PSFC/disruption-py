function year = create_year_column;
db = set_database('d3drdb');

shot = cell2mat(fetch(db, ['select shot from disruption_warning ' ...
  'order by shot, time']));

shotmin = min(shot);

result = fetch(db, ['select shot, time_of_shot from summaries ' ...
  'where shot >= ' num2str(shotmin) ' order by shot']);
shotlist = int32(cell2mat(result(:,1)));

year = NaN(length(shot), 1);
for i=1:length(shotlist);
  time_of_shot_string = char(result(i,2));
  yearval = str2num(time_of_shot_string(1:4));
  indices = find(ismember(shot, shotlist(i)));
  if length(indices) > 0;
    year(indices) = yearval;
  end;
end;
