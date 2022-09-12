%% Get Data
variables = {'li', 'shot', 'time', 'time_until_disrupt', 'ip', ...
  'ip_error', 'dipprog_dt', 'intentional_disruption', ...
  'other_hardware_failure', 'power_supply_railed'};

% Turn off annoying warning messages that tend to be generated when using
% the SQL java driver

warning_status = warning;
warning('off')

db_retrieve = set_database_d3d('d3drdb');

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

shotlist_2015 = dlmread(['d3d_2015_shotlist.txt']);

% Get the indices of time slices from the flattop of non-disruptive shots
% in the 2015 shotlist

define_indices;

indices_flattop_disrupt_2015 = intersect(...
    indices_flattop_disrupt_in_flattop, find(ismember(shot, ...
    shotlist_2015)));

shotlist_disrupt_2015 = unique(shot(indices_flattop_disrupt_2015));
nshots = length(shotlist_disrupt_2015);

%% Plot Figure

count_disruptions_li_rise = 0;
shots_with_li_rise_during_long_flattop = NaN(nshots,1);

figure('menubar', 'none', 'toolbar', 'none', 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','Position',[0.1,0.53,0.80,0.37]);
clf;
plot(NaN(1),NaN(1),'.'); xlim([-1,0]); ylim([0.8,1.6]);
set(gca,'fontsize',12);
xlabel('Time Until Disrupt [s]','fontsize',14);
ylabel('l_i','fontsize',14);
title('l_i behavior during flattop of disruptive shots in 2015', ...
  'fontsize',15);
hold on;

for ishot = 1:nshots;
  fprintf(1,'Processing shot %6i  (%4i/%4i  %6.2f%%)\n', ...
    shotlist_disrupt_2015(ishot), ishot, nshots, ishot/nshots*100);

  td = find(time_until_disrupt==0 & shot==shotlist_disrupt_2015(ishot));
  
  %---------------------------------------------------------------------%
  % Get early and late li stats for current shot

  late_li = li(find(time(td)-time < 0.1 & ...
      shot==shotlist_disrupt_2015(ishot)));
  early_li = li(find(time(td)-time > 0.8 & time(td)-time < 1.0 & ...
      shot==shotlist_disrupt_2015(ishot)));

  early_li_mean = mean(early_li,'omitnan');
  early_li_std = std(early_li,'omitnan');
  late_li_mean = mean(late_li,'omitnan');
  late_li_std = std(late_li,'omitnan');
  %---------------------------------------------------------------------%
  
  shot_flattop_indices = intersect(...
      find(shot==shotlist_disrupt_2015(ishot)), ...
      indices_flattop_disrupt_in_flattop);

  if (early_li_mean < 1.0 && early_li_std < 0.05 && ...
      late_li_mean > 1.1 && late_li_mean < 1.3 && ...
      late_li_std < 0.1); % shows behavior
    count_disruptions_li_rise = count_disruptions_li_rise + 1;
    plot(-time_until_disrupt(shot_flattop_indices), li(shot_flattop_indices), '.-', ...
        'color', [255, 69, 0]/256);
  else
    plot(-time_until_disrupt(shot_flattop_indices), li(shot_flattop_indices), '.-', ...
        'color', [1,1,1]*0.70);
  end;
end;

fprintf(1, 'Count of disruptions with desired behavior of li = %i\n', ...
  count_disruptions_li_rise);

hold off;
print('li_behavior_non_disrupt.png', '-dpng');
