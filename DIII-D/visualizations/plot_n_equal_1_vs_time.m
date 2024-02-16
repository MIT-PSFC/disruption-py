variables = {'n_equal_1_normalized', 'time_until_disrupt', 'shot', ...
  'time', 'ip', 'ip_error', 'dipprog_dt', 'intentional_disruption', ...
  'other_hardware_failure', 'power_supply_railed'};

% Turn off annoying warning messages that tend to be generated when using
% the SQL java driver

warning_status = warning;
warning('off')

db_retrieve = set_database('d3drdb');

for ivar = 1:length(variables);
  varname = char(variables(ivar));
  if exist(varname, 'var') == 0;
    charlength = length(varname);
    fprintf(1, 'retrieving %s\n', varname);
    result_retrieve = fetch(db_retrieve, ['select ' varname ...
      ' from disruption_warning order by shot, time']);
    eval([varname '= cell2mat(result_retrieve);']);
  end;
end;
close(db_retrieve);

define_indices;
jj = indices_flattop_disrupt_in_flattop;

% Get the actual DIII-D 2015 shotlist used in our PPCF paper

shotlist_2015 = dlmread(['/fusion/projects/disruption_warning/' ...
  'paper_programs/PPCF2018/d3d_2015_shotlist.txt']);

% Get the indices of time slices from disruptive shots in the 2015 shotlist

jjj = intersect(jj, find(ismember(shot, shotlist_2015)));

% Get the list of disruptive shots in the 2015 shotlist

shotlist_2015_disrupt = unique(shot(jjj));

% Plot

set(gcf, 'units', 'centimeters', 'position', [5, 18, 40, 10]);
clf;
plot([-1,0],[NaN,NaN]);
xlim([-1,0]); ylim([0,1.e-3]);
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [s]', 'fontsize', 14);
ylabel('B_p^{n=1}/B_{tor}', 'fontsize', 14);
title('n=1 magnitude vs time until disrupt; DIII-D 2015', 'fontsize', 15)
hold on;
for ishot = 1:length(shotlist_2015_disrupt);
  jjjj = intersect(jjj, find(shot == shotlist_2015_disrupt(ishot)));
  plot(-time_until_disrupt(jjjj), n_equal_1_normalized(jjjj), '.-');
  pause(0.1);
end;
hold off;

% Return warning status back to its original setting

warning(warning_status);
