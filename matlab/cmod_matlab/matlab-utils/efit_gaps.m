% Checks for time periods in shots from DIII-D database where 
% multiple EFIT variables are all missing, or NaN.
%
% Author: Kevin Montes
% Date: Dec 2019
%___________________________________________________________________

%% Grab EFIT data, shot/time indices, and list of shots
% Check for SQL data, and get it if not present
if ~exist('vars')==1
    disp('Loading data ...')
    SQL_retrieve
else
    disp('Data already loaded ...')
end

% Get a list of all shots missing EFITs during some portion of flattop
flattop_indx = union(indices_flattop_disrupt_in_flattop, ...
    indices_flattop_no_disrupt);
missing = find(isnan(vars.q95) & isnan(vars.li) & isnan(vars.beta_p));
missing_in_flattop = intersect(flattop_indx,missing);
shots_with_nans = unique(shot(missing_in_flattop));

% Use this section ONLY IF excluding intentional disruptions & hardware
% failures from shots_with_nans
fail_indices = union(union(find(vars.intentional_disruption==1), ...
    find(vars.other_hardware_failure==1)),find(vars.power_supply_railed==1));
failed_shots = unique(shot(fail_indices));
shots_with_nans = setdiff(shots_with_nans,failed_shots);

% Divide shots with gaps into 'disruption' and 'non-disruption' categories
disruptions = unique(shot(indices_flattop_disrupt_in_flattop));
non_disrupts = unique(shot(indices_flattop_no_disrupt));
disrupts_with_nans = intersect(shots_with_nans,disruptions);
non_disrupts_with_nans = intersect(shots_with_nans,non_disrupts);
included_disruptions = setdiff(disruptions,failed_shots);

%% Loop through shots with gaps to find gap lengths

num_gaps = nan(size(shots_with_nans));
largest_gap_size = nan(size(shots_with_nans));
smallest_gap_size = nan(size(shots_with_nans));
first_gap_tud = nan(size(disrupts_with_nans));
all_missing = logical(zeros(size(shots_with_nans)));
with_end_gap = logical(zeros(size(shots_with_nans)));
last_gaptime = nan(size(shots_with_nans));
norm_time = time;
N = length(shots_with_nans);
for i=1:length(shots_with_nans)
    disp(['Doing shot ' num2str(i) ' of ' num2str(N)])
    shot_num = shots_with_nans(i); % shot number
    shot_indx = intersect(flattop_indx,find(shot==shot_num));
    nan_indx = intersect(shot_indx,missing); % find indices of missing data
    if length(shot_indx) == length(nan_indx)
        all_missing(i) = 1;
        continue
    end
    [first_pt,last_pt] = locate_gaps(nan_indx); % see local function
    gaptimes = diff(time(nan_indx([first_pt; last_pt]))); % in [s]
    if ~isempty(gaptimes)
        if nan_indx(last_pt(end)) == shot_indx(end)
            with_end_gap(i) = 1;
            last_gaptime(i) = diff(time(nan_indx([first_pt(end); last_pt(end)])));
        end
        num_gaps(i) = length(gaptimes); 
        smallest_gap_size(i) = min(gaptimes);
        largest_gap_size(i) = max(gaptimes);
        if ismember(shot_num,disrupts_with_nans)
            jj = find(disrupts_with_nans==shot_num);
            first_gap_tud(jj) = time_until_disrupt(nan_indx(first_pt(1)));
        end
    end
    t = time(shot_indx);
    norm_time(shot_indx) = (t-t(1))/range(t);
end

%% Display info graphically

% How many shots don't have any EFITs at all?
num_missing = length(find(all_missing));
num_shots = length(shots_with_nans);
disp([num2str(num_missing) ' of ' num2str(num_shots) ...
    ' shots with gaps (' num2str(round(num_missing/num_shots*100,2)) ...
    '%) are missing EFITs for all time slices'])

only_some_missing = find(ismember(shot,shots_with_nans(~all_missing)));
indx_of_interest = intersect(missing_in_flattop,only_some_missing);
disr_interest = intersect(indx_of_interest,indices_disrupt);
nondisr_interest = intersect(indx_of_interest,indices_no_disrupt);

% What does the missing time slice distribution look like for each shot? 
figure();
plot(-time_until_disrupt(indx_of_interest),shot(indx_of_interest),'.')
xlabel('Time Until Disrupt [s]')
ylabel('Shot #')
title('Time Records Missing EFITs')
ylim([min(shot),max(shot)])

figure(); hold on;
plot(norm_time(nondisr_interest),shot(nondisr_interest),'.')
plot(norm_time(disr_interest),shot(disr_interest),'r.')
xlabel('$(t-t_0)/(t_f-t_0)$','Interpreter','latex')
ylabel('Shot #')
title('Time Records Missing EFITs')
legend('Non-disrupt','Disrupted')
ylim([min(shot),max(shot)])

% Histogram of time_until_disrupt for gap start time
figure();
histogram(-first_gap_tud,[-1:0.005:0],'normalization','probability')
xlabel('Time Until Disrupt of First EFIT Gap [s]')
ylabel('Fraction of flattop shots')

% Histogram of number of EFIT gaps
figure();
disr_gap = zeros(size(shots_with_nans));
disr_gap(find(ismember(shots_with_nans,disrupts_with_nans))) = 1;
histogram(num_gaps(~all_missing),'normalization','probability')
xlabel('No. of EFIT gaps')
ylabel('Fraction of flattop shots')

% Histogram of smallest gap time
figure();
histogram(smallest_gap_size(~all_missing),'normalization','probability')
xlabel('Smallest Gap [s]')
ylabel('Fraction of flattop shots')

% Histogram of largest gap time
figure();
histogram(largest_gap_size(~all_missing),'normalization','probability')
xlabel('Largest Gap [s]')
ylabel('Fraction of flattop shots')

%% Plot cumulative histogram of gap time for disruptions
% only disruptive shots with a gap at end
jj = intersect(find(ismember(shots_with_nans,disrupts_with_nans)), ...
    find(ismember(shots_with_nans,shots_with_nans(with_end_gap)))); 
figure(); hold on;
h = histogram(-[last_gaptime(jj); nan(length(included_disruptions) - ...
    length(jj),1)]*1e3, -fliplr([0:0.01:max(last_gaptime)]*1e3), ...
    'normalization','cdf')
thresh = h.Values(max(find(h.BinEdges<-50)));
plot(xlim,[thresh thresh],'k--')
xlabel('Last gap time [ms]')
ylabel('Fraction of disruptions')
xlim([-1,0]*1e3)
title('Cumulative Distribution of EFIT Gap Times Near Disruption')

%% Test

indx = [10];
[first_pt,last_pt] = locate_gaps(indx);
gaptimes = diff(time(nan_indx([first_pt; last_pt]))); % in [s]

%% Local functions

function [first_pt,last_pt] = locate_gaps(indx)
% Given indices 'indx' of records with missing data, this function 
% identifies end points of gaps in the signal

if iscolumn(indx); 
    indx = indx'; 
elseif ~isrow(indx); 
    disp('Error... Input indices must be a vector');
    return
end

t = diff(indx) == 1;
y = [t,false];
x = xor(y,[false,t]);
ii = cumsum(~(x|y) + y.*x);
jumpIndx = find(diff(ii));
first_pt = [1 jumpIndx+1];
last_pt = [jumpIndx numel(ii)];
to_keep = find(first_pt~=last_pt);
first_pt = first_pt(to_keep);
last_pt = last_pt(to_keep);

end
