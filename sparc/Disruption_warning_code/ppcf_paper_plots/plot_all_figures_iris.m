%-------------------------------------------------------------%
% This routine generates each of the plots in the file
% 'summary_of_possible_figures.pdf' using a standard font size
% for each title and axis label, as specified at the beginning 
% of the code.
%-------------------------------------------------------------%

addpath('/fusion/projects/disruption_warning/matlab_programs')

%% Set fontsizes

title_size = 16;
axis_label_size = 14;

%% Start DIII-D section
%-----------------------------------------------------------%
% Connect to DIII-D disruption warning database (may need 
% to modify)
%-----------------------------------------------------------%
dbconn = set_database('d3drdb');

%% Fetch needed data from DIII-D database

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

%% Select subset of blessed shots from 2015 DIII-D campaign

shotlist = int32(dlmread('d3d_2015_shotlist.txt'));
indices_d3d_2015 = find(ismember(shot,shotlist));

dipprog_dt = dipprog_dt(indices_d3d_2015);
intentional_disruption = intentional_disruption(indices_d3d_2015);
ip = ip(indices_d3d_2015);
ip_error = ip_error(indices_d3d_2015);
other_hardware_failure = other_hardware_failure(indices_d3d_2015);
power_supply_railed = power_supply_railed(indices_d3d_2015);
shot = shot(indices_d3d_2015);
li = li(indices_d3d_2015);
time = time(indices_d3d_2015);
time_until_disrupt = time_until_disrupt(indices_d3d_2015);

clearvars dbconn

define_indices;

%% Start li histogram for DIII-D

ii = indices_flattop_no_disrupt;
jj = indices_flattop_disrupt_in_flattop;
jj_far = intersect(jj, find(time_until_disrupt > 0.35));
jj_near = intersect(jj, find(time_until_disrupt <= 0.35 & ...
  time_until_disrupt > 0.015)); % 2 ms 'black window'

bins = 0.00 : 0.02 : 2.00;

figure;
histogram(li(ii), bins, 'normalization', 'probability', ...
  'facecolor', 'blue');
xlim([0.5, 2.0]);
ylim([0.00, 0.12]);
set(gca, 'fontsize', 12);
xlabel('l_i', 'fontsize', axis_label_size);
ylabel('probability histogram', 'fontsize', axis_label_size);
title('Histograms of l_i for DIII-D', 'fontsize', title_size);
hold on;
histogram(li(jj_far), bins, 'normalization', 'probability', ...
  'facecolor', 'green');
histogram(li(jj_near), bins, 'normalization', 'probability', ...
  'facecolor', 'red');
legend('non-disruptive', 'time until disrupt > 350 ms', ...
  'time until disrupt < 350 ms', 'location', 'northeast');
hold off;

%% Start li vs. time_until_disrupt for DIII-D

ii = indices_flattop;
figure; hold on;
for i=1:length(shotlist)
    indices = intersect(ii,find(shot==shotlist(i)));
    plot(-time_until_disrupt(indices),li(indices),'b.-')
end
xlabel('Time until disrupt [s]','fontsize', axis_label_size)
ylabel('l_i','fontsize', axis_label_size)
title('l_i for all flattop disruption times (DIII-D)', ...
    'fontsize', title_size)
xlim([-1,0])
ylim([0.8,1.6])
hold off;

%% Get mean li data for DIII-D

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

figure;
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
close(gcf); % Close the generated histogram (don't need it)

clearvars -except mean_li_d3d xaxis_d3d axis_label_size ...
    title_size

%% Start C-Mod Section
%-----------------------------------------------------------%
% Connect to C-Mod disruption warning database (may need
% need to modify)
%-----------------------------------------------------------%

dbconn = set_database_cmod('logbook');

%% Fetch needed data from C-Mod database

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

%% Select subset of blessed shots from 2015 C-Mod campaign

shotlist = int32(dlmread('cmod_2015_shotlist.txt'));
indices_cmod_2015 = find(ismember(shot,shotlist));

dipprog_dt = dipprog_dt(indices_cmod_2015);
ip = ip(indices_cmod_2015);
ip_error = ip_error(indices_cmod_2015);
power_supply_railed = power_supply_railed(indices_cmod_2015);
shot = shot(indices_cmod_2015);
li = li(indices_cmod_2015);
time = time(indices_cmod_2015);
time_until_disrupt = time_until_disrupt(indices_cmod_2015);

clearvars dbconn

define_indices;

%% Produce histograms of internal inductance (li) for C-Mod

ii = indices_flattop_no_disrupt;
jj = indices_flattop_disrupt_in_flattop;
jj_far = intersect(jj, find(time_until_disrupt > 0.040));
jj_near = intersect(jj, find(time_until_disrupt <= 0.040 & ...
  time_until_disrupt > 0.002)); % 2 ms 'black window'

bins = 0.00 : 0.02 : 2.00;
figure;
histogram(li(ii), bins, 'normalization', 'probability', ...
  'facecolor', 'blue');
xlim([0.5, 2.0]);
ylim([0.00, 0.11]);
set(gca, 'fontsize', 12);
xlabel('l_i', 'fontsize', axis_label_size);
ylabel('probability histogram', 'fontsize', axis_label_size);
title('Histograms of l_i for Alcator C-Mod', 'fontsize', ...
    title_size);
hold on;
histogram(li(jj_far), bins, 'normalization', 'probability', ...
  'facecolor', 'green');
histogram(li(jj_near), bins, 'normalization', 'probability', ...
  'facecolor', 'red');
legend('non-disruptive', 'time until disrupt > 40 ms', ...
  'time until disrupt < 40 ms', 'location', 'northwest');
hold off;

%% Plot li for all flattop disruption times

ii = indices_flattop;

figure; hold on;
for i=1:length(shotlist)
    indices = intersect(ii,find(shot==shotlist(i)));
    plot(-time_until_disrupt(indices),li(indices),'b.-')
end
xlabel('Time until disrupt [s]','fontsize', axis_label_size)
ylabel('l_i','fontsize', axis_label_size)
title('l_i for all flattop disruption times (C-Mod)', ...
    'fontsize', title_size)
xlim([-1,0])
ylim([0.8,1.6])
hold off;

%% Get mean li data for C-Mod

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
close(gcf)

clearvars -except axis_label_size title_size ...
    mean_li_cmod mean_li_d3d xaxis_cmod xaxis_d3d

%% Generate mean li plot and zoom inset for C-Mod & DIII-D

figure;
plot(xaxis_cmod, mean_li_cmod, 'sb');
hold on;
plot(xaxis_d3d, mean_li_d3d, 'sr');
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [s]', 'fontsize', axis_label_size);
ylabel('mean(l_i)', 'fontsize', axis_label_size);
title('Mean of l_i distribution vs time before disrupt', ...
  'fontsize', title_size);
legend('C-Mod','DIII-D','location','northwest')
hold off;

figure;
plot(xaxis_cmod(9:end),mean_li_cmod(9:end),'sb')
xlabel('Time until disrupt [s]','fontsize',axis_label_size)
ylabel('mean(l_i)','fontsize', axis_label_size)
title('Mean of l_i distribution vs time before disrupt (C-Mod)', ...
    'fontsize',title_size)