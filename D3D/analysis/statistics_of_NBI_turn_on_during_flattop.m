dbconn = set_database('d3drdb');  % Connect to SQL database

% Get variables that are need for the routine "define_indices.m":

if exist('shot','var') == 0;
  result = fetch(dbconn,['select shot from ' ...
    'disruption_warning order by shot, time']);
  shot = int32(cell2mat(result));
end;
if exist('time','var') == 0;
  result = fetch(dbconn,['select time from ' ...
    'disruption_warning order by shot, time']);
  time = cell2mat(result);
end;
if exist('time_until_disrupt','var') == 0;
  result = fetch(dbconn,['select time_until_disrupt from ' ...
    'disruption_warning order by shot, time']);
  time_until_disrupt = cell2mat(result);
end;
if exist('ip','var') == 0;
  result = fetch(dbconn,['select ip from ' ...
    'disruption_warning order by shot, time']);
  ip = cell2mat(result);
end;
if exist('ip_error','var') == 0;
  result = fetch(dbconn,['select ip_error from ' ...
    'disruption_warning order by shot, time']);
  ip_error = cell2mat(result);
end;
if exist('dipprog_dt','var')==0;
  result = fetch(dbconn,['select dipprog_dt from ' ...
    'disruption_warning order by shot, time']);
  dipprog_dt = cell2mat(result);
end;
if exist('power_supply_railed','var')==0;
  result = fetch(dbconn,['select power_supply_railed from ' ...
    'disruption_warning order by shot, time']);
  power_supply_railed = cell2mat(result);
end;
if exist('intentional_disruption','var')==0;
  result = fetch(dbconn,['select intentional_disruption from ' ...
    'disruption_warning order by shot, time']);
  intentional_disruption = cell2mat(result);
end;
if exist('other_hardware_failure','var')==0;
  result = fetch(dbconn,['select other_hardware_failure from ' ...
    'disruption_warning order by shot, time']);
  other_hardware_failure = cell2mat(result);
end;

% Get the NBI power:

if exist('p_nbi','var')==0;
  result = fetch(dbconn,['select p_nbi from ' ...
    'disruption_warning order by shot, time']);
  p_nbi = cell2mat(result);
end;

close(dbconn);

define_indices;  % needed to select flattop times

indices_nbi_on = find(p_nbi >= 0.5e6);
shotlist = unique(shot);
nshots = length(shotlist);
flattop_start_time = NaN(1,nshots);
nbi_turn_on_time = NaN(1,nshots);
nbi_flattop_start_delay = NaN(1, nshots);

for ishot = 1:nshots;
  indices_this_shot = find(shot == shotlist(ishot));
  indices_flattop_this_shot = ...
    intersect(indices_flattop, indices_this_shot);
  indices_nbi_on_this_shot = ...
    intersect(indices_nbi_on, indices_this_shot);
  if length(indices_flattop_this_shot) > 0;
    flattop_start_time(ishot) = time(indices_flattop_this_shot(1));
  end;
  if length(indices_nbi_on_this_shot) > 0;
    nbi_on_time = time(indices_nbi_on_this_shot(1));
  end;
  nbi_turn_on_time(ishot) = nbi_on_time;
  nbi_flattop_start_delay(ishot) = nbi_on_time - flattop_start_time(ishot);
end;

save('nbi_flattop_start_delay.mat', 'nbi_flattop_start_delay', ...
  'nbi_turn_on_time', 'flattop_start_time', 'shotlist');

figure;
histogram(flattop_start_time, [-2 : 0.1 : 10], ...
  'normalization', 'probability');
xlim([0, 5]);
set(gca, 'fontsize', 12);
xlabel('Start time of flattop [s]', 'fontsize', 14);
ylabel('probability histogram', 'fontsize', 14);
title('Statistics of NBI turn-on time during flattop', 'fontsize', 15)

print('flattop_start_time.png', '-dpng');

figure;
histogram(nbi_turn_on_time, [-2 : 0.1 : 10], ...
  'normalization', 'probability');
xlim([-2.5, +10.5]);
set(gca, 'fontsize', 12);
xlabel('Start time of NBI [s]', 'fontsize', 14);
ylabel('probability histogram', 'fontsize', 14);
title('Statistics of NBI turn-on time during flattop', 'fontsize', 15)

print('nbi_flattop_start_time.png', '-dpng');

figure;
histogram(nbi_flattop_start_delay, [-2 : 0.1 : 10], ...
  'normalization', 'probability');
xlim([-2.5, +3.5]);
set(gca, 'fontsize', 12);
xlabel('Time between start of flattop and start of NBI [s]', ...
  'fontsize', 14);
ylabel('probability histogram', 'fontsize', 14);
title('Statistics of NBI turn-on time during flattop', 'fontsize', 15)

print('nbi_flattop_start_delay.png', '-dpng');
