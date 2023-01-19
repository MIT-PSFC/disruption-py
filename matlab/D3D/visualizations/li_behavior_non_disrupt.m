% Get necessary data.  (Many of these are needed for the "define_indices"
% routine)

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

indices_flattop_no_disrupt_2015 = intersect(...
    indices_flattop_no_disrupt, find(ismember(shot, shotlist_2015)));

shotlist_no_disrupt_2015 = unique(shot(indices_flattop_no_disrupt_2015));
nshots = length(shotlist_no_disrupt_2015);

% Plot Figure

figure('menubar', 'none', 'toolbar', 'none', 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','Position',[0.1,0.53,0.80,0.37]);
clf;
plot(NaN, NaN, '.'); xlim([0,8]); ylim([0.8,1.6]);
set(gca,'fontsize',12);
%xlabel('Adjusted time [s]','fontsize', 14);
xlabel('Time [s]','fontsize', 14);
ylabel('l_i','fontsize',14);
title('l_i behavior during flattop of non-disruptive shots in 2015', ...
  'fontsize',15);
hold on;

for loopnum = 1:2;  % I need to do this so that the orange lines on the
                    % graph don't get covered over by the more numerous
                    % gray lines.

count_no_disrupt_with_long_enough_flattop = 0;
count_no_disrupt_li_rise_during_long_enough_flattop = 0;
shots_with_li_rise_during_long_enough_flattop = NaN(nshots,1);

for ishot = 1:nshots;
  fprintf(1,'Processing shot %6i  (%4i/%4i  %6.2f%%)\n', ...
    shotlist_no_disrupt_2015(ishot), ishot, nshots, ishot/nshots*100);

  flattop_indices = intersect(indices_flattop_no_disrupt, ...
    find(shot == shotlist_no_disrupt_2015(ishot)));
  flattop_indices = intersect(flattop_indices, ...
    find(time > 1.4)); % J(r) is still peaking up during first 1.4 s of shot
  flattop_times = time(flattop_indices);
  li_values = li(flattop_indices);

  flattop_duration = max(flattop_times) - min(flattop_times);
  if flattop_duration < 1.0; continue; end;

  count_no_disrupt_with_long_enough_flattop = ...
    count_no_disrupt_with_long_enough_flattop + 1;

% Use a FOR loop to sweep a one-second wide window along the flattop, and
% examine each one-second wide window to determine if the desired li
% temporal behavior occurs in the window.  (For each shot, start the sweep
% at the end of flattop, and move backwards in 0.5 s steps.)
  
  li_behavior = 0;
  for endtime = max(flattop_times) : -0.500 : min(flattop_times)+1.0;
    indices_window = find((flattop_times >= endtime - 1.0) & ...
      (flattop_times <= endtime));
    indices_early = find((flattop_times >= endtime - 1.0) & ...
      (flattop_times <= endtime - 0.8));
    indices_late = find((flattop_times >= endtime - 0.1) & ...
      (flattop_times <= endtime));
    
% Get early and late li stats for current shot

    early_li = li_values(indices_early);
    late_li = li_values(indices_late);

    early_li_mean = mean(early_li,'omitnan');
    early_li_std = std(early_li,'omitnan');
    late_li_mean = mean(late_li,'omitnan');
    late_li_std = std(late_li,'omitnan');
    
    if (early_li_mean < 1.0 && early_li_std < 0.05 && ...
        late_li_mean > 1.1 && late_li_mean < 1.3 && ...
        late_li_std < 0.1); % shows behavior
      li_behavior = 1;
%     plot(flattop_times(indices_window) - endtime, ...
%       li_values(indices_window), '.-', ...
%       'color', [255, 69, 0]/256);
    else
%     plot(flattop_times(indices_window) - endtime, ...
%       li_values(indices_window), '.-', ...
%       'color', [1,1,1]*0.70);
    end
  end;
  if (li_behavior == 1); % Only count each shot once
    count_no_disrupt_li_rise_during_long_enough_flattop = ...
      count_no_disrupt_li_rise_during_long_enough_flattop + 1;
    if loopnum == 2;
      plot(flattop_times, li_values, 'color', [255, 69, 0]/256);
    end;
  else;
    if loopnum == 1;
      plot(flattop_times, li_values, 'color', [1,1,1]*0.70);
    end;
  end;
end;
end;

fprintf(1, 'Count of non-disruptions with any length flattop = %i\n', nshots);
fprintf(1, 'Count of non-disruptions with long enough flattop = %i\n', ...
  count_no_disrupt_with_long_enough_flattop);
fprintf(1, 'Count of non-disruptions with desired behavior of li = %i\n', ...
  count_no_disrupt_li_rise_during_long_enough_flattop);

hold off;
print('li_behavior_non_disrupt.png', '-dpng');
