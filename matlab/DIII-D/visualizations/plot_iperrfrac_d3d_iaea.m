addpath('/fusion/projects/disruption_warning/matlab_programs/');

dbconn = set_database('d3drdb');

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

% Get the fractional Ip error

if exist('ip_error','var')==0;
  result = fetch(dbconn,['select ip_error from ' ...
    'disruption_warning order by shot, time']);
  ip_error = cell2mat(result);
end;
if exist('ip_prog','var')==0;
  result = fetch(dbconn,['select ip_prog from ' ...
    'disruption_warning order by shot, time']);
  ip_prog = cell2mat(result);
end;

close(dbconn);

ip_error_frac = ip_error ./ ip_prog;

define_indices;  % needed to select flattop times
ii = indices_flattop_no_disrupt;
jj = indices_flattop_disrupt_in_flattop;

figure('units', 'centimeters', 'position', [10, 10, 15.0, 6.0]);
plot(-smear(time_until_disrupt(jj),.002), ip_error_frac(jj), '.b');
xlim([-0.1, 0]);
ylim([-5,0]);
set(gca, 'fontsize', 15);
xlabel('Time until disrupt [s]', 'fontsize', 14);
ylabel('(I_p error)/I_p', 'fontsize', 14);
text(0.02, 0.90, 'DIII-D', 'fontsize', 16, 'fontweight', 'bold', ...
  'units', 'normalized', 'backgroundcolor', 'w', 'margin', 1);
print('iperrfrac_vs_time_until_disrupt_d3d_iaea.png', '-dpng');
