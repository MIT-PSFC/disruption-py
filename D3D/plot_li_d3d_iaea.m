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

% Get the internal inductance signal:

if exist('li','var')==0;
  result = fetch(dbconn,['select li from ' ...
    'disruption_warning order by shot, time']);
  li = cell2mat(result);
end;

close(dbconn);

define_indices;  % needed to select flattop times

ii = indices_flattop_no_disrupt;
jj = indices_flattop_disrupt_in_flattop;

% Randomly reduce the sampling rate of near-disruption time slices (every 2
% ms) to match the sampling rate of far-from-disrupt time slices (every 25
% ms).  This is solely for aesthetic reasons.

[indices_thinned] = thin_points(time_until_disrupt(jj), 0.1001, 2/25);
jj = jj(indices_thinned);

% There are so many points on the DIII-D plot that we can't see the
% important behavior.  There, randomly reduce the number of all points
% within 1.0 s of a disruption by a factor of 5 (because we're only
% plotting the period within 1.0 s before disruption).

[indices_thinned] = thin_points(time_until_disrupt(jj), 1.001, 0.2);
jj = jj(indices_thinned);

figure('units', 'centimeters', 'position', [10, 10, 15.0, 6.0]);
plot(-smear(time_until_disrupt(jj),.002), li(jj), '.b');
xlim([-1.0, 0]);
ylim([0.75, 1.4]);
set(gca, 'fontsize', 15);
xlabel('Time until disrupt [s]', 'fontsize', 14);
ylabel('l_i', 'fontsize', 14);
text(0.02, 0.90, 'DIII-D', 'fontsize', 16, 'fontweight', 'bold', ...
  'units', 'normalized', 'backgroundcolor', 'w', 'margin', 1);
print('li_vs_time_until_disrupt_d3d_iaea.png', '-dpng');
