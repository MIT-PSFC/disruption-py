dbconn = set_database_cmod('logbook');

% Fetch needed data from DIII-D database

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
if exist('li','var')==0;
  result = fetch(dbconn,['select li from ' ...
    'disruption_warning order by shot, time']);
  li = cell2mat(result);
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
close(dbconn);

define_indices;

indices_2015 = find(shot>1150000000 & shot<1160000000);

ii = intersect(indices_flattop_no_disrupt,indices_2015);
jj = intersect(indices_flattop_disrupt_in_flattop,indices_2015);
jj_far = intersect(jj, find(time_until_disrupt > 0.040));
jj_near = intersect(jj, find(time_until_disrupt <= 0.040 & ...
  time_until_disrupt > 0.002)); % 2 ms 'black window'

bins = 0.00 : 0.02 : 2.00;

histogram(li(ii), bins, 'normalization', 'probability', ...
  'facecolor', 'blue');
xlim([0.5, 2.0]);
ylim([0.00, 0.11]);
set(gca, 'fontsize', 12);
xlabel('l_i', 'fontsize', 14);
ylabel('probability histogram', 'fontsize', 14);
title('Histograms of l_i for Alcator C-Mod', 'fontsize', 16);
hold on;
histogram(li(jj_far), bins, 'normalization', 'probability', ...
  'facecolor', 'green');
histogram(li(jj_near), bins, 'normalization', 'probability', ...
  'facecolor', 'red');
legend('non-disruptive', 'time until disrupt > 40 ms', ...
  'time until disrupt < 40 ms', 'location', 'northwest');
hold off;