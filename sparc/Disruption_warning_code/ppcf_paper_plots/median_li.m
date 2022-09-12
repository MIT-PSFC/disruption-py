% Median li analysis

addpath('/home/granetz/JRT_2016/disruption_warning_database');

%% Set fontsizes

%title_size = 16; Chose to exclude titles from plots
axis_label_size = 14;

%% Start DIII-D section
%-----------------------------------------------------------%
% Connect to DIII-D disruption warning database (may need 
% to modify)
%-----------------------------------------------------------%
dbconn = set_database_d3d('d3drdb');

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
if exist('Te_HWHM','var')==0;
  result = fetch(dbconn,['select Te_HWHM from ' ...
    'disruption_warning order by shot, time']);
  Te_HWHM = cell2mat(result);
end;
if exist('Greenwald_fraction','var')==0;
  result = fetch(dbconn,['select Greenwald_fraction from ' ...
    'disruption_warning order by shot, time']);
  Greenwald_fraction = cell2mat(result);
end;
if exist('n_equal_1_normalized','var')==0;
  result = fetch(dbconn,['select n_equal_1_normalized' ...
    ' from disruption_warning order by shot, time']);
  n_equal_1_normalized = cell2mat(result);
end;
if exist('beta_p','var')==0;
  result = fetch(dbconn,['select beta_p' ...
    ' from disruption_warning order by shot, time']);
  beta_p = cell2mat(result);
end;
if exist('q95','var')==0;
  result = fetch(dbconn,['select q95' ...
    ' from disruption_warning order by shot, time']);
  q95 = cell2mat(result);
end;
if exist('v_loop','var')==0;
  result = fetch(dbconn,['select v_loop' ...
    ' from disruption_warning order by shot, time']);
  v_loop = cell2mat(result);
end;
if exist('Wmhd','var')==0;
  result = fetch(dbconn,['select Wmhd' ...
    ' from disruption_warning order by shot, time']);
  Wmhd = cell2mat(result);
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
time = time(indices_d3d_2015);
time_until_disrupt = time_until_disrupt(indices_d3d_2015);

ip_error_frac = ip_error./ip;
beta_p = beta_p(indices_d3d_2015);
Greenwald_fraction = Greenwald_fraction(indices_d3d_2015);
q95 = q95(indices_d3d_2015);
li = li(indices_d3d_2015);
Te_HWHM = Te_HWHM(indices_d3d_2015);
n_equal_1_normalized = n_equal_1_normalized(indices_d3d_2015);
v_loop = v_loop(indices_d3d_2015);
Wmhd = Wmhd(indices_d3d_2015);

clearvars dbconn

define_indices;

%% Start median li plots

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
ii_00 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0 & time_until_disrupt <= 0.1));

xaxis_d3d = [-1.05, -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35,...
    -0.25, -0.15,-0.05];

li_median = NaN(1,12);

li_median(1) = median(li(ii_far),'omitnan');
li_median(2) = median(li(ii_10),'omitnan');
li_median(3) = median(li(ii_09),'omitnan');
li_median(4) = median(li(ii_08),'omitnan');
li_median(5) = median(li(ii_07),'omitnan');
li_median(6) = median(li(ii_06),'omitnan');
li_median(7) = median(li(ii_05),'omitnan');
li_median(8) = median(li(ii_04),'omitnan');
li_median(9) = median(li(ii_03),'omitnan');
li_median(10) = median(li(ii_02),'omitnan');
li_median(11) = median(li(ii_01),'omitnan');
li_median(12) = median(li(ii_00),'omitnan');

median_li_d3d = li_median(2:end);

%% Clear workspace for next section

clearvars -except median_li_d3d xaxis_d3d axis_label_size ...
    title_size

%% Start C-Mod Section
%-----------------------------------------------------------%
% Connect to C-Mod disruption warning database (may need
% need to modify)
%-----------------------------------------------------------%

dbconn = set_database('logbook');

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

close(dbconn);

%% Select subset of blessed shots from 2015 C-Mod campaign

shotlist = int32(dlmread('cmod_2015_shotlist.txt'));
indices_cmod_2015 = find(ismember(shot,shotlist));

dipprog_dt = dipprog_dt(indices_cmod_2015);
ip = ip(indices_cmod_2015);
shot = shot(indices_cmod_2015);
li = li(indices_cmod_2015);
time = time(indices_cmod_2015);
time_until_disrupt = time_until_disrupt(indices_cmod_2015);

clearvars dbconn

define_indices;

%% Start the mean li plot for C-Mod

ii_far = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 1.1));
ii_18 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 1 & time_until_disrupt <= 1.1));
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

xaxis = [-1.05,-0.95,-0.85,-0.75,-0.65,-0.55,-0.45,-0.35,-0.25, ...
    -0.19, -0.17, -0.15, -0.13, -0.11, -0.09, -0.07, -0.05, -0.03]% -0.01];

li_median = NaN(1,19);

li_median(1) = median(li(ii_far),'omitnan');
li_median(2) = median(li(ii_18),'omitnan');
li_median(3) = median(li(ii_17),'omitnan');
li_median(4) = median(li(ii_16),'omitnan');
li_median(5) = median(li(ii_15),'omitnan');
li_median(6) = median(li(ii_14),'omitnan');
li_median(7) = median(li(ii_13),'omitnan');
li_median(8) = median(li(ii_12),'omitnan');
li_median(9) = median(li(ii_11),'omitnan');
li_median(10) = median(li(ii_10),'omitnan');
li_median(11) = median(li(ii_09),'omitnan');
li_median(12) = median(li(ii_08),'omitnan');
li_median(13) = median(li(ii_07),'omitnan');
li_median(14) = median(li(ii_06),'omitnan');
li_median(15) = median(li(ii_05),'omitnan');
li_median(16) = median(li(ii_04),'omitnan');
li_median(17) = median(li(ii_03),'omitnan');
li_median(18) = median(li(ii_02),'omitnan');
li_median(19) = median(li(ii_01),'omitnan');

median_li_cmod = li_median(2:end);
xaxis_cmod = xaxis;
close(gcf)

clearvars -except axis_label_size title_size ...
    median_li_cmod median_li_d3d xaxis_cmod xaxis_d3d

figure;
plot(xaxis_cmod, median_li_cmod, 'sb');
hold on;
plot(xaxis_d3d, median_li_d3d, 'sr');
set(gca, 'fontsize', 12);
xlabel('time\_until\_disrupt [s]', 'fontsize', axis_label_size, ...
    'FontName','MathJax_Typewriter');
ylabel('median(li)', 'fontsize', axis_label_size,...
    'FontName','MathJax_Typewriter');
%title('Mean of l_i distribution vs time before disrupt', ...
%  'fontsize', title_size);
legend('C-Mod','DIII-D','location','northwest')
hold off;

figure;
plot(xaxis_cmod(10:end),median_li_cmod(10:end),'sb')
xlabel('time\_until\_disrupt [s]','fontsize',axis_label_size, ...
    'FontName','MathJax_Typewriter');
ylabel('median(li)','fontsize', axis_label_size, ...
    'FontName','MathJax_Typewriter');
%title('Mean of l_i distribution vs time before disrupt (C-Mod)', ...
%    'fontsize',title_size)