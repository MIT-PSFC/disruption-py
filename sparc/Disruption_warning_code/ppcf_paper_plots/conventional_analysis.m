%----------------------------------------------------------------------%
%
% This script implements a 'conventional analysis' of the 2015 DIII-D 
% disruption warning database using thresholds on the li parameters. It 
% is meant to respond to the peer review critique for the PPCF submission
% 
% Author: Kevin Montes
% Date: 04/02/2018
%
%----------------------------------------------------------------------%

addpath('/home/granetz/JRT_2016/disruption_warning_database');
dbconn = set_database_d3d('d3drdb');
axis_label_size = 14;
title_size = 16;

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
if exist('dipprog_dt','var')==0;
  result = fetch(dbconn,['select dipprog_dt from ' ...
    'disruption_warning order by shot, time']);
  dipprog_dt = cell2mat(result);
end;
if exist('li','var')==0;
  result = fetch(dbconn,['select li from ' ...
    'disruption_warning order by shot, time']);
  li = cell2mat(result);
end;

%% Select subset of blessed shots from 2015 DIII-D campaign

shotlist = int32(dlmread('d3d_2015_shotlist.txt'));
indices_d3d_2015 = find(ismember(shot,shotlist));

dipprog_dt = dipprog_dt(indices_d3d_2015);
ip = ip(indices_d3d_2015);
shot = shot(indices_d3d_2015);
time = time(indices_d3d_2015);
time_until_disrupt = time_until_disrupt(indices_d3d_2015);
li = li(indices_d3d_2015);

clearvars dbconn

define_indices;
ii = indices_flattop_no_disrupt;
jj = indices_flattop_disrupt_in_flattop;

%% Histograms for li

time_thresholds = [0.2:0.02:0.6]; % in seconds
overlap = NaN(size(time_thresholds));
for t = 1:length(time_thresholds)
    thresh = time_thresholds(t);
    jj_far = intersect(jj, find(time_until_disrupt > thresh));
    jj_near = intersect(jj, find(time_until_disrupt <= thresh & ...
      time_until_disrupt > 0.015)); % 2 ms 'black window'
    jj_stable = union(ii,jj_far);

    bins = 0.00 : 0.02 : 2.00;
    figure;
    h_far = histogram(li(jj_far), bins, 'normalization', ...
        'probability');
    xlim([0.5, 2.0]);
    ylim([0.00, 0.12]);
    set(gca, 'fontsize', 12);
    xlabel('l_i', 'fontsize', axis_label_size);
    ylabel('probability histogram', 'fontsize', axis_label_size);
    %title('Histograms of l_i for DIII-D', 'fontsize', title_size);
    hold on;
    h_near = histogram(li(jj_near), bins, 'normalization', 'probability');
    legend(['Far from disruption (>' num2str(1000*thresh) ' ms)'], ...
      ['Near disruption (<' num2str(1000*thresh) ' ms'], ...
      'location', 'northeast');
    hold off;

    overlap(t) = sum(min([h_near.Values;h_far.Values])); % record integral
end
[overlap_val,min_indx] = min(overlap);
disp(['Threshold w/ smallest histogram overlap: ' ...
    num2str(time_thresholds(min_indx))])

%% New Section


%% Start li vs time plot for DIII-D (non-disruptive shots)

figure;
plot(time(ii),li(ii),'.-')
ylim([0,4])

%% Start li vs. time_until_disrupt for DIII-D (disruptive shots)

disrupted_shotlist = unique(shot(find(~isnan(time_until_disrupt))));
nColors = 8;
cmap = linspace(0,1,nColors+1);
cmap = [cmap;cmap;cmap]';

bins = linspace(0.8,1.6,nColors);

figure; hold on;
red_list = [];
for i=1:length(disrupted_shotlist)
    indices = intersect(jj,find(shot==disrupted_shotlist(i)));
    tTemp = time_until_disrupt(indices);
    liTemp = li(indices);
    [~,tindx] = min(abs(tTemp-1));
    bool = liTemp(tindx) < bins;
    [~,indx] = min(abs(bool-1));
    if liTemp(tindx) > 0.9 && liTemp(tindx) < 1.0
        map = [0.8500 0.3250 0.0980];
        [~,tindx] = min(abs(tTemp-0.02));
        if liTemp(tindx)>1.1 & liTemp(tindx)<1.3
            red_list = [red_list,disrupted_shotlist(i)];
        end
    else
        map = cmap(indx,:);
    end
    plot(-tTemp,liTemp,'.-', ...
        'Color', map)% [0 0.4470 0.7410])
end
xlabel('time\_until\_disrupt [s]','fontsize', axis_label_size,...
    'FontName','MathJax_Typewriter')
ylabel('li','fontsize', axis_label_size,'FontName','MathJax_Typewriter')
title('l_i for all flattop disruption times (DIII-D)', ...
    'fontsize', title_size)
xlim([-1,0])
ylim([0.8,1.6])
set(gcf,'units','centimeters','position',[1,1,20,10]);
box on
hold off;

figure; hold on;
map = [0.8500 0.3250 0.0980];
for i=1:length(red_list)
    indices = intersect(jj,find(shot==red_list(i)));
    plot(-time_until_disrupt(indices),li(indices),'.-','Color',map)
end
xlim([-1,0])
ylim([0.8,1.6])
set(gcf,'units','centimeters','position',[1,1,20,10]);
box on
hold off;
disp([num2str(length(red_list)) '/' ...
    num2str(length(disrupted_shotlist)) ' shots highlighted'])

%% Test time and li thresholds for classification

% Define thresholds to be sweeped over
time_thresholds = 0.35; % in seconds
li_thresholds = [0.9:0.01:1.5];

% Get flattop data and initiate class vectors
td = time_until_disrupt(indices_flattop);
t = time(indices_flattop);
s = shot(indices_flattop);
X = li(indices_flattop);
y = zeros(size(indices_flattop));
y_pred = zeros(size(indices_flattop));
FN = NaN(length(time_thresholds),length(li_thresholds));
FP = NaN(size(FN));
TN = NaN(size(FN));
TP = NaN(size(FN));

% Sweep over each time and li threshold, redefining classes and making
% predictions for each before recording each corresponding F1 score
for i=1:length(time_thresholds)
    y = zeros(size(indices_flattop));
    y(find(td < time_thresholds(i))) = 1; % define classes
    for j=1:length(li_thresholds)
        y_pred = zeros(size(indices_flattop));
        y_pred(find(X > li_thresholds(j))) = 1; % predict
        TP(i,j) = length(find(y == y_pred & y==1)); % true positives
        FN(i,j) = length(find(y ~= y_pred & y==1)); % false negatives
        FP(i,j) = length(find(y ~= y_pred & y==0)); % false positives
        TN(i,j) = length(find(y == y_pred & y==0)); % true negatives
    end
end

% Calculate F1 score for each threshold pair
F1_score = 2*TP./(2*TP+FP+FN);
recall = TP./(TP+FN);
precision = TP./(TP+FP);
figure;
pcolor(time_thresholds,li_thresholds,F1_score')
colorbar
xlabel('Time threshold [s]'); ylabel('li threshold');
title('Recall Using Univariate li Analysis')
%}