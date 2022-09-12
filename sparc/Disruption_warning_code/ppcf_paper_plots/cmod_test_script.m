%% Start C-Mod Section
%-----------------------------------------------------------%
% Connect to C-Mod disruption warning database (may need
% need to modify)
%-----------------------------------------------------------%

dbconn = set_database('logbook');
axis_label_size = 14;

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
if exist('n_equal_1_normalized','var')==0;
  result = fetch(dbconn,['select n_equal_1_normalized from ' ...
    'disruption_warning order by shot, time']);
  n_equal_1_normalized = cell2mat(result);
end

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
n_equal_1_normalized = n_equal_1_normalized(indices_cmod_2015);

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
ylim([0.00, 0.12]);
set(gca, 'fontsize', 12);
xlabel('l_i', 'fontsize', axis_label_size);
ylabel('probability histogram', 'fontsize', axis_label_size);
%title('Histograms of l_i for Alcator C-Mod', 'fontsize', ...
%    title_size);
hold on;
histogram(li(jj_far), bins, 'normalization', 'probability', ...
  'facecolor', 'green');
histogram(li(jj_near), bins, 'normalization', 'probability', ...
  'facecolor', 'red');
legend('non-disruptive', 'time until disrupt > 40 ms', ...
  'time until disrupt < 40 ms', 'location', 'northwest');
hold off;

%% Histograms for n_equal_1_mode on C-Mod

% Histogram for n_equal_1_normalized

bins = 0:3e-5:3e-3;
figure;
histogram(n_equal_1_normalized(ii), bins, 'normalization', ...
  'probability', 'facecolor', 'blue');
%xlim([0, 7e-4]);
%ylim([0.00, 0.09]);
set(gca, 'fontsize', 12);
xlabel('n\_equal\_1\_normalized', 'fontsize', axis_label_size);
ylabel('probability histogram', 'fontsize', axis_label_size);
%title('Histograms of l_i for C-Mod', 'fontsize', title_size);
hold on;
histogram(n_equal_1_normalized(jj_far), bins, ...
  'normalization', 'probability', 'facecolor', 'green');
histogram(n_equal_1_normalized(jj_near), bins, ...
  'normalization', 'probability', 'facecolor', 'red');
legend('non-disruptive', 'time until disrupt > 40 ms', ...
  'time until disrupt < 40 ms', 'location', 'northeast');
hold off;

%% Plot li for all flattop disruption times on C-Mod

ii = indices_flattop;
nColors = 6;
cmap = linspace(0,1,nColors+1);
cmap = [cmap;cmap;cmap]';

bins = linspace(1.1,1.6,nColors);

figure; hold on;
for i=1:length(shotlist)
    indices = intersect(ii,find(shot==shotlist(i)));
    bool = li(indices(1)) < bins;
    [~,indx] = min(abs(bool-1));
    plot(-time_until_disrupt(indices),li(indices),'.-', ...
        'Color',cmap(indx,:));
end
xlabel('Time until disrupt [s]','fontsize', axis_label_size)
ylabel('l_i','fontsize', axis_label_size)
%title('l_i for all flattop disruption times (C-Mod)', ...
%    'fontsize', title_size)
xlim([-1,0])
ylim([0.8,1.6])
set(gcf,'units','centimeters','position',[1,1,20,10]);
box on
hold off;