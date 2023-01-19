variables = {'li', 'shot', 'time', 'time_until_disrupt', 'ip', ...
  'ip_error', 'dipprog_dt', 'intentional_disruption', ...
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
clearvars db_retrieve result_retrieve varname charlength variables;

% Get the actual DIII-D 2015 shotlist used in our PPCF paper

shotlist_2015 = dlmread(['/fusion/projects/disruption_warning/' ...
  'paper_programs/PPCF2018/d3d_2015_shotlist.txt']);

% Get the indices of time slices from the flattop of non-disruptive shots
% in the 2015 shotlist

define_indices;

indices_flattop_no_disrupt_2015 = intersect(indices_flattop_no_disrupt, ...
  find(ismember(shot, shotlist_2015)));

shotlist_no_disrupt_2015 = unique(shot(indices_flattop_no_disrupt_2015));
nshots = length(shotlist_no_disrupt_2015);

count_no_disrupt_with_flattop = 0;
count_no_disrupt_with_long_flattop = 0;
count_no_disrupt_li_rise_during_long_flattop = 0;
shots_with_li_rise_during_long_flattop = NaN(nshots,1);

figure('menubar', 'none', 'toolbar', 'none', 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','InnerPosition',[0.1,0.53,0.80,0.37]);
clf;
plot(NaN(1),NaN(1),'.'); xlim([0,7]); ylim([0.8,1.6]);
set(gca,'fontsize',12);
xlabel('Time [s]','fontsize',14);
ylabel('l_i','fontsize',14);
title('l_i behavior during flattop of non-disruptive shots in 2015', ...
  'fontsize',15);
hold on;

for ishot = 1:nshots;
  fprintf(1,'Processing shot %6i  (%4i/%4i  %6.2f%%)\n', ...
    shotlist_no_disrupt_2015(ishot), ishot, nshots, ishot/nshots*100);

  count_no_disrupt_with_flattop = count_no_disrupt_with_flattop + 1;

  flattop_indices = intersect(find(shot == ...
    shotlist_no_disrupt_2015(ishot)), indices_flattop_no_disrupt);
  flattop_times = time(flattop_indices);
  flattop_duration = flattop_times(end) - flattop_times(1);

  if flattop_duration < 1.25; continue; end;

  count_no_disrupt_with_long_flattop = count_no_disrupt_with_long_flattop + 1;

  flattop_li = li(flattop_indices);

  plot(flattop_times, flattop_li, '.-', 'color', [1,1,1]*0.70);

  indices_early_flattop = find(flattop_times <= flattop_times(1) + 1.0);
  indices_late_flattop = find(flattop_times >= flattop_times(end) - 0.25);

  early_flattop_li = flattop_li(indices_early_flattop);
  late_flattop_li = flattop_li(indices_late_flattop);

  early_flattop_li_mean = mean(early_flattop_li);
  early_flattop_li_std = std(early_flattop_li);
  late_flattop_li_mean = mean(late_flattop_li);
  late_flattop_li_std = std(late_flattop_li);

  if (early_flattop_li_mean < 1.0 && early_flattop_li_std < 0.05 && ...
      late_flattop_li_mean > 1.1 && late_flattop_li_mean < 1.3 && ...
      late_flattop_li_std < 0.1);
    count_no_disrupt_li_rise_during_long_flattop = ...
      count_no_disrupt_li_rise_during_long_flattop + 1;
    shots_with_li_rise_during_long_flattop( ...
      count_no_disrupt_li_rise_during_long_flattop) = ...
      shotlist_no_disrupt_2015(ishot);
  end;
end;

shots_with_li_rise_during_long_flattop = ...
  shots_with_li_rise_during_long_flattop(find(~isnan( ...
  shots_with_li_rise_during_long_flattop)));

fprintf(1, 'Count of non-disruptions with flattop = %i\n', ...
  count_no_disrupt_with_flattop);
fprintf(1, 'Count of non-disruptions with long enough flattop = %i\n', ...
  count_no_disrupt_with_long_flattop);
fprintf(1, 'Count of non-disruptions with desired behavior of li = %i\n', ...
  count_no_disrupt_li_rise_during_long_flattop);

for ishot = 1:length(shots_with_li_rise_during_long_flattop);
  qq = find(shot == shots_with_li_rise_during_long_flattop(ishot));
  flattop_indices = intersect(indices_flattop_no_disrupt, qq);
  plot(time(flattop_indices), li(flattop_indices), '.-', 'color', ...
    [255, 69, 0]/256);
end;

hold off;
print('li_behavior_non_disrupt.png', '-dpng');
