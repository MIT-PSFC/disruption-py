addpath(['/fusion/projects/disruption_warning' ...
    '/matlab_programs']);

%-----------------------------------------------------------%
% Connect to DIII-D disruption warning database (may need 
% to modify)
%-----------------------------------------------------------%
dbconn = set_database_d3d('d3drdb');

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

% select subset of blessed 392 discharges

shotlist = int32(dlmread('shotlist_rea_blessed.txt'));
%
shotlist=dlmread(['/fusion/projects/disruption_warning/matlab_programs/' ...
  'shotlist_rea_blessed.txt']);
shotlist = int32(shotlist);
%
subset_blessed = find(ismember(shot,shotlist));

dipprog_dt = dipprog_dt(subset_blessed);
intentional_disruption = intentional_disruption(subset_blessed);
ip = ip(subset_blessed);
ip_error = ip_error(subset_blessed);
other_hardware_failure = other_hardware_failure(subset_blessed);
power_supply_railed = power_supply_railed(subset_blessed);
shot = shot(subset_blessed);
li = li(subset_blessed);
time = time(subset_blessed);
time_until_disrupt = time_until_disrupt(subset_blessed);

clearvars dbconn shotlist

define_indices;

% Start histogram for DIII-D

ii_far = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 1.1));
ii_10 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 1.0 & time_until_disrupt <= 1.1));
ii_09 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.9 & time_until_disrupt <= 1.0));
ii_08 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.8 & time_until_disrupt <= 0.9));
ii_07 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.7 & time_until_disrupt <= 0.8));
ii_06 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.6 & time_until_disrupt <= 0.7));
ii_05 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.5 & time_until_disrupt <= 0.6));
ii_04 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.4 & time_until_disrupt <= 0.5));
ii_03 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.3 & time_until_disrupt <= 0.4));
ii_02 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.2 & time_until_disrupt <= 0.3));
ii_01 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.1 & time_until_disrupt <= 0.2));

xaxis_d3d = [-1.05, -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15];

bins = 0 : .02 : 2;
bin_centers = bins(1:end-1) + diff(bins)/2;

h0 = histogram(li(ii_far), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
set(gca, 'fontsize', 12);
xlabel('l_i', 'fontsize', 14);
ylabel('probability', 'fontsize', 14);
title('Histograms of l_i (DIII-D) at various times before disrupt', 'fontsize', 16);
hold on;

h10 = histogram(li(ii_10), bins, 'normalization', 'probability' , ...
  'facecolor', 'black');
h9 = histogram(li(ii_09), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');
h8 = histogram(li(ii_08), bins, 'normalization', 'probability' , ...
  'facecolor', 'red');
h7 = histogram(li(ii_07), bins, 'normalization', 'probability' , ...
  'facecolor', 'yellow');
h6 = histogram(li(ii_06), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
h5 = histogram(li(ii_05), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');
h4 = histogram(li(ii_04), bins, 'normalization', 'probability' , ...
  'facecolor', 'red');
h3 = histogram(li(ii_03), bins, 'normalization', 'probability' , ...
  'facecolor', 'yellow');
h2 = histogram(li(ii_02), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
h1 = histogram(li(ii_01), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');

% Re-plot li(ii_far) so it can be seen on top of all the others.

histogram(li(ii_far), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');

hold off;

mean_li = NaN(1,11);

mean_li(1) = xcentroid(bin_centers, h0.Values);
mean_li(2) = xcentroid(bin_centers, h10.Values);
mean_li(3) = xcentroid(bin_centers, h9.Values);
mean_li(4) = xcentroid(bin_centers, h8.Values);
mean_li(5) = xcentroid(bin_centers, h7.Values);
mean_li(6) = xcentroid(bin_centers, h6.Values);
mean_li(7) = xcentroid(bin_centers, h5.Values);
mean_li(8) = xcentroid(bin_centers, h4.Values);
mean_li(9) = xcentroid(bin_centers, h3.Values);
mean_li(10) = xcentroid(bin_centers, h2.Values);
mean_li(11) = xcentroid(bin_centers, h1.Values);

mean_li_d3d = mean_li(2:end);
clearvars -except mean_li_d3d xaxis_d3d

%-----------------------------------------------------------%
% Connect to C-Mod disruption warning database (may need
% need to modify)
%-----------------------------------------------------------%

dbconn = set_database_cmod('logbook');

% Fetch needed data from C-Mod database

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

close(dbconn);

define_indices;

%subselect shots from 2015 campaign
indices_2015 = find(shot>1150000000 & shot<1160000000);
indices_flattop_disrupt_in_flattop = intersect(...
    indices_flattop_disrupt_in_flattop,indices_2015);

% Start the histogram for C-Mod

ii_far = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 1));
ii_17 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.9 & time_until_disrupt <= 1));
ii_16 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.8 & time_until_disrupt <= 0.9));
ii_15 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.7 & time_until_disrupt <= 0.8));
ii_14 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.6 & time_until_disrupt <= 0.7));
ii_13 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.5 & time_until_disrupt <= 0.6));
ii_12 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.4 & time_until_disrupt <= 0.5));
ii_11 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.3 & time_until_disrupt <= 0.4));
ii_10 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.2 & time_until_disrupt <= 0.3));
ii_09 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.18 & time_until_disrupt <= 0.2));
ii_08 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.16 & time_until_disrupt <= 0.18));
ii_07 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.14 & time_until_disrupt <= 0.16));
ii_06 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.12 & time_until_disrupt <= 0.14));
ii_05 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.10 & time_until_disrupt <= 0.12));
ii_04 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.08 & time_until_disrupt <= 0.10));
ii_03 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.06 & time_until_disrupt <= 0.08));
ii_02 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.04 & time_until_disrupt <= 0.06));
ii_01 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.02 & time_until_disrupt <= 0.04));

xaxis = [-0.95,-0.85,-0.75,-0.65,-0.55,-0.45,-0.35,-0.25, ...
    -0.19, -0.17, -0.15, -0.13, -0.11, -0.09, -0.07, -0.05, -0.03]% -0.01];

bins = 1 : .005 : 1.8;
bin_centers = bins(1:end-1) + diff(bins)/2;

figure;
h0 = histogram(li(ii_far), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
set(gca, 'fontsize', 12);
xlabel('l_i', 'fontsize', 14);
ylabel('probability', 'fontsize', 14);
title('Histograms of l_i (C-Mod) at various times before disrupt', 'fontsize', 16);
hold on;

h17 = histogram(li(ii_17), bins, 'normalization', 'probability' , ...
  'facecolor', 'yellow');
h16 = histogram(li(ii_16), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
h15 = histogram(li(ii_15), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');
h14 = histogram(li(ii_14), bins, 'normalization', 'probability' , ...
  'facecolor', 'red');
h13 = histogram(li(ii_13), bins, 'normalization', 'probability' , ...
  'facecolor', 'yellow');
h12 = histogram(li(ii_12), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
h11 = histogram(li(ii_11), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');
h10 = histogram(li(ii_10), bins, 'normalization', 'probability' , ...
  'facecolor', 'black');
h9 = histogram(li(ii_09), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');
h8 = histogram(li(ii_08), bins, 'normalization', 'probability' , ...
  'facecolor', 'red');
h7 = histogram(li(ii_07), bins, 'normalization', 'probability' , ...
  'facecolor', 'yellow');
h6 = histogram(li(ii_06), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
h5 = histogram(li(ii_05), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');
h4 = histogram(li(ii_04), bins, 'normalization', 'probability' , ...
  'facecolor', 'red');
h3 = histogram(li(ii_03), bins, 'normalization', 'probability' , ...
  'facecolor', 'yellow');
h2 = histogram(li(ii_02), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
h1 = histogram(li(ii_01), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');

% Re-plot li(ii_far) so it can be seen on top of all the others.

histogram(li(ii_far), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');

hold off;

mean_li = NaN(1,18);

mean_li(1) = xcentroid(bin_centers, h0.Values);
mean_li(2) = xcentroid(bin_centers, h17.Values);
mean_li(3) = xcentroid(bin_centers, h16.Values);
mean_li(4) = xcentroid(bin_centers, h15.Values);
mean_li(5) = xcentroid(bin_centers, h14.Values);
mean_li(6) = xcentroid(bin_centers, h13.Values);
mean_li(7) = xcentroid(bin_centers, h12.Values);
mean_li(8) = xcentroid(bin_centers, h11.Values);
mean_li(9) = xcentroid(bin_centers, h10.Values);
mean_li(10) = xcentroid(bin_centers, h9.Values);
mean_li(11) = xcentroid(bin_centers, h8.Values);
mean_li(12) = xcentroid(bin_centers, h7.Values);
mean_li(13) = xcentroid(bin_centers, h6.Values);
mean_li(14) = xcentroid(bin_centers, h5.Values);
mean_li(15) = xcentroid(bin_centers, h4.Values);
mean_li(16) = xcentroid(bin_centers, h3.Values);
mean_li(17) = xcentroid(bin_centers, h2.Values);
mean_li(18) = xcentroid(bin_centers, h1.Values);

mean_li_cmod = mean_li(2:end);
xaxis_cmod = xaxis;

clearvars -except mean_li_cmod mean_li_d3d xaxis_cmod xaxis_d3d

figure;
plot(xaxis_cmod, mean_li_cmod, 'sb');
hold on;
plot(xaxis_d3d, mean_li_d3d, 'sr');
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [s]', 'fontsize', 14);
ylabel('mean(l_i)', 'fontsize', 14);
title('Mean of l_i distribution vs time before disrupt', ...
  'fontsize', 16);
legend('C-Mod','DIII-D','location','northwest')

figure;
plot(xaxis_cmod(9:end),mean_li_cmod(9:end),'sb')
xlabel('Time until disrupt [s]','fontsize',14)
ylabel('mean(l_i)','fontsize',14)
title('Mean of l_i distribution vs time before disrupt (C-Mod)', ...
    'fontsize',16)