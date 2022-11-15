function shotlist = create_list_of_plasmas_in_2015;
db = set_database('d3drdb');

% Specifications for minimum pulse_length and minimum ip were decided in a
% conversation I had with Nick Eidietis on 2016/02/25.

result = fetch(db, ['select shot, time_of_shot from summaries ' ...
  'where pulse_length > 0.1 and abs(ip) > 400e3 ' ...
  'order by shot']);
shot = int32(cell2mat(result(:,1)));
time_of_shot = cell2mat(result(:,2));
indices = NaN(length(shot), 1);
for i=1:length(shot);
  time_of_shot_string = char(time_of_shot(i));
  indices(i) = strcmpi(time_of_shot_string(1:4), '2015');
end;
indices = find(indices);
shotlist = shot(indices);
