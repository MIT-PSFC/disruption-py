% The program "get_n1rms_d3d.m", written by Jinxiang Zhu, returns values
% for n1rms_normalized in units of gauss/tesla, and these values were
% loaded into our disruption_warning database in early October 2019.  I
% prefer to have the values be dimensionless.  So I wrote this little
% program to multiply all the values in the "n1rms_normalized" column by
% 10^-4.
%
% RSG  2019/10/07

fprintf(1, 'Fetching dbkey, n1rms_normalized from database\n');

db = set_database('d3drdb');
result = fetch(db, ['select dbkey, n1rms_normalized from ' ...
  'disruption_warning order by dbkey']);
record_keys = cell2mat(result(:,1));
n1rms_normalized = cell2mat(result(:,2));
  
fprintf(1, 'Done fetching dbkey, n1rms_normalized from database\n');

for i = 1:length(record_keys);
  if (isnan(n1rms_normalized(i)) || ...
      isinf(n1rms_normalized(i)));
    value = {''};
  else;
    value = {n1rms_normalized(i) * 1e-4};
  end;
  update(db, 'disruption_warning', {'n1rms_normalized'}, value, ...
    ['where dbkey = ' num2str(record_keys(i),'%i')]);
end;
