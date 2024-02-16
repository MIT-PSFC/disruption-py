function [list_correct_FD, list_early_FD, list_miss_FD, list_late_FD, ...
    list_correct_ND, list_fp_ND, ...
    list_pd_d, list_pd_nd, list_pd_fp, list_pd_m, list_pd_l] = ...
    threshold_algorithm_v2(shotlist, parameters, threshold_matrix, t_lower, t_upper)

% I'm updating this to plot a threshold scan as in threshold_scan.m and its
% several versions.


% Here is my attempt at a disruption warning algorithm using thresholds
% I've found while exploring the trends in data as different parameters
% approach a disruption. 
%
% Inputs:
%   shotlist: shots to be tested
%   parameters: string (cell) array of N parameters
%   threshold matrix: 2xN matrix, each row corresponds to one parameter and
%       the 2 columns to the lower and upper thresholds, respectively
%   t_lower = lower bound for correctly ID disruption (usually 10ms)
%             so afterwards it's too late to mitigate
%   t_upper = upper bound for correctly ID disruption (usually 20ms)
%             Also separates (pre)disruptive data from the "non" disruptive
%             data of a disruption
%
% Alex Tinguely 160411

% The following code was stolen from plot_disruption_data.m.
db = set_database('logbook');

shotlist = unique(sort(shotlist)); % put the shots in order from smallest to largest
                                   % and remove redundant shots if necessary

parameters = lower(strtrim(parameters)); % parameter to be fetched

% Get data from ALL shots in the shotlist (provided they are in the db)
data = cell2mat(fetch(db, ['select shot, time, time_until_disrupt, ip, ', ... 
    'ip_error, dipprog_dt, ', strjoin(parameters, ','), ' from disruption_warning ',...
    'where shot between ', num2str(shotlist(1)), ' and ', num2str(shotlist(end)),...
    ' order by shot, time']));
% Note the following parameters:
%   data(:,1) = shot
%   data(:,2) = time
%   data(:,3) = time_until_disrupt
%   data(:,4) = ip
%   data(:,5) = ip_error = ip - ip_prog
%   data(:,6) = dipprog_dt
%   data(:,6) = parameter

% Now we only want data from shots in the shotlist so I check to see if the
% shots as listed in column 1 of data (redundant) are included in our
% shotlist. If a shot is, it gets a 1 (for all timeslices); if not, we get
% a 0. This is effective in only using shots from the shotlist that are
% also in the database.

is_data_in_shotlist = ismember(data(:,1), shotlist); % boolean array
data = data(is_data_in_shotlist, :);

% We are only interested in data during the plasma current flattop, where
% Ip_prog > 100kA (should actually be around 1MA), and dIpprog_dt < 60kA/s
% (should actually be 0, but normal ramp up/down rates are >100kA/s.)

ip_prog = data(:,4) - data(:,5); % = ip - ip_error 
is_flattop_1 = abs(data(:,6)) < 6e4; % dipprog_dt < 60kA/s (
is_flattop_2 = abs(ip_prog) > 1e5; % ip_prog > 100kA
is_flattop = and(is_flattop_1, is_flattop_2);

data = data(is_flattop,:); % just data during the flattop for both non/disruptions

% Now we want to split up our data into disruptions and non-disruptions.
% First, we check if the time_until_disrupt is NaN, which is indicative of
% a non-disruption. 

t_is_nan = isnan(data(:,3)); % time_until_disrupt is not NaN

nd_data = data(t_is_nan,:); % non-disruptions
d_data = data(~t_is_nan,:); % disruptions

% Find the unique shots in each list.
nd_shots = sort(unique(nd_data(:,1))); % non-disruption shotlist
d_shots = sort(unique(d_data(:,1))); % disruption shotlist

% We also want to be careful of "non-disruptions" that last for less
% (approximately) the full two seconds. We have now restricted our data to
% those timeslices during the flattop portion (~1-1.5s). Let's make sure our
% *non-disruptions* don't end before 1s. (This could be upped.)
% NOTE: These shots are rare, so it probably gets rid of 1/1000.

for i = 1:length(nd_shots)
    
    bool = ismember(nd_data(:,1),nd_shots(i)); % find timeslices for ith nd_shot
    t = nd_data(bool,2); % gets the time array for the ith nd_shot
    
    if max(t) < 1 % if the nd_shot time is less than 1s
        nd_shots(i) = -1; % set the shot to -1 to identify after the loop
    end
end

% Now remove all of the -1's from nd_shots
nd_shots(nd_shots == -1) = [];

% Now we want only the data in nd_data corresponding to nd_shots.

bool2 = ismember(nd_data(:,1),nd_shots);
nd_data = nd_data(bool2,:);

% We only want disruptions that actually occur during the flattop.
% Currently we have flattop data, but there is the possibility that
% time_until_disrupt = 0 occurs after the flattop. We will label flattop
% disruptions with "fd" or "FD". The other disruptions will be
% labeled as "rd" or "RD", because in this case, these other disruptions occur in
% the Rampdown. (All rampup disruptions were cut because they have no
% flattop data.)

% boolean if time_until_disrupt == 0 so disruption has occurred in flattop
is_t_until_disrupt_0 = d_data(:,3) == 0;

% find shot numbers where disrupton occurs during flattop
shots_t0 = unique(d_data(is_t_until_disrupt_0, 1));
% get a boolean array of the timeslices for these shots
is_shot_with_t0 = ismember(d_data(:,1),shots_t0);
FD_data = d_data(is_shot_with_t0,:);

% The Rampdown disruptions don't have a time_until_disrupt == 0. However,
% they still may have times within t_d (usually 20ms) of a disruption.

RD_data = d_data(~is_shot_with_t0,:);

% Just renaming non-disruption data.
ND_data = nd_data;

% Also get rid of disruption data with t_until_disrupt < 0.

FD_data(FD_data(:,3) < 0,:) = [];

% Get the shots in each
ND_shots = sort(unique(ND_data(:,1)));
RD_shots = sort(unique(RD_data(:,1)));
FD_shots = sort(unique(FD_data(:,1)));

% Now we will walk through each shot and count the number of false
% positives, early warnings, etc.

%------------------Flattop Disruptions FD---------------%

% initialize the number number of definitions for *flattop disruptions* FD
n_early_FD = 0; % early warnings
n_correct_FD = 0; % correctly identified disruptions within t_lower, t_upper
n_late_FD = 0; % late warning
n_miss_FD = 0; % missed disruption

% Create lists for shots within the above categories.
list_correct_FD = [];
list_early_FD = [];
list_miss_FD = [];
list_late_FD = [];

% First consider the disruption shots.
for i = 1:length(FD_shots)
    
    % get boolean array for all timeslices of the ith flattop disruption shot.
    % This uses FD_data from above which inludes all timeslices for the
    % disruption (even t_until_disrupt > 20ms)
    bool_FD = ismember(FD_data(:,1), FD_shots(i));
    
    temp_data = FD_data(bool_FD, :); % all data for just the ith shot
    
    % initialize boolean of zeros
    label_FD = boolean(zeros(size(temp_data,1),1)); 
    
    % For each of the parameters (starting at column 7), we see if the data
    % is past the thresholds as provided in threshold_matrix. Then we use
    % the or() function to compare with the other parameters. This
    % "label_FD" then indicates whether the datapoint was considered to be a
    % disruption (1 = yes, 0 = no).
    
    for j = 1:length(parameters)
        
        is_past_threshold = or(temp_data(:,6+j) <= threshold_matrix(j,1),...
            temp_data(:,6+j) >= threshold_matrix(j,2));
        
        label_FD = or(label_FD, is_past_threshold);
        
    end
    
    % times before disruption at which a disruption warning is given
    t = temp_data(label_FD, 3); 
    
    if isempty(t) % if no disruption warning was given, it was missed!
        n_miss_FD = n_miss_FD + 1; 
        list_miss_FD = [list_miss_FD, int32(FD_shots(i))];
    else      
        tmax = max(t);
        
        % if the disruption warning is given before t_d s before the disruption
        if tmax > t_upper 
            % we consider it too early (during normal plasma)
            n_early_FD = n_early_FD +1; 
            list_early_FD = [list_early_FD, int32(FD_shots(i))];
            
        % if the disruption warning is given between t_lower and t_upper
        elseif tmax >= t_lower 
            % it is correctly identified in time
            n_correct_FD = n_correct_FD + 1;
            list_correct_FD = [list_correct_FD, int32(FD_shots(i))];
         
        % if the warning is within t_lower of the disruption, it's too late    
        else 
            n_late_FD = n_late_FD + 1;
            list_late_FD = [list_late_FD, int32(FD_shots(i))];
            
        end        
    end
end

%------------------Non-Disruptions ND--------------------%

% For non-disruptions ND, we only have
n_fp_ND = 0; % false positive
n_correct_ND = 0; % correctly identified non-disruptions

% So create lsits
list_correct_ND = [];
list_fp_ND = [];

% Next look at the non-disruption shots
for i = 1:length(ND_shots)
       
    % This uses ND_data from above, but still have to select out the
    % data for each test shot.
    bool_ND2 = ismember(ND_data(:,1), ND_shots(i));
    temp_data = ND_data(bool_ND2, :); % all data for just the ith shot
    
     % initialize boolean of zeros
    label_ND = boolean(zeros(size(temp_data,1),1)); 
    
    for j = 1:length(parameters)
        
        is_past_threshold = or(temp_data(:,6+j) <= threshold_matrix(j,1),...
            temp_data(:,6+j) >= threshold_matrix(j,2));
        
        label_ND = or(label_ND, is_past_threshold);
        
    end
    
    % if any disruption warning is given, it's a false positive
    if sum(label_ND) > 0 
        n_fp_ND = n_fp_ND + 1;
        list_fp_ND = [list_fp_ND, int32(ND_shots(i))];

    % or else it is a correctly idenitified non-disruption    
    else
        n_correct_ND = n_correct_ND + 1; 
        list_correct_ND = [list_correct_ND, int32(ND_shots(i))];
    end
end

%--------------------Rampdown Disruptions RD-------------%

% Now we look at shots that are rampdown disruptions, where they disrupt
% *after* the flattop. Thus, for the most part, they are just like
% non-disruptions. However, if the disruption occurs relatively soon after
% the flattop end (within t_upper to t_lower s), then we can sort of treat this like a
% disruption. If this is the case, can we even predict these disruptions?

% First let's separate the RD_data into those shots that have
% time_until_disrupt data that is within t_upper s and those that do not.
% The latter are basically non-disruptions. The former we are interested to
% see if they could still be predicted even though the disruption is not in
% the flattop.

contains_t_less_than_t_upper = RD_data(:,3) <= t_upper;

% This data is "disruption-like". (probably very few shots)
RD_D_data = RD_data(contains_t_less_than_t_upper, :);
RD_D_shots = sort(unique(RD_D_data(:,1)));

% The rest of the data is "non-disruption-like".
RD_ND_shots = sort(unique(setdiff(RD_shots, RD_D_shots)));
bool5 = ismember(RD_data(:,1), RD_ND_shots);
RD_ND_data = RD_data(bool5,:);

% Now we follow a similar procedure for the ND-like shots.

n_fp_RD_ND = 0; % false positive
n_correct_RD_ND = 0; % correctly identified non-disruptions

% So create lsits
list_correct_RD_ND = [];
list_fp_RD_ND = [];

% Look at the rampdown non-disruptive shots
for i = 1:length(RD_ND_shots)
       
    bool_RD_ND = ismember(RD_ND_data(:,1), RD_ND_shots(i));
    temp_data = RD_ND_data(bool_RD_ND, :); % all data for just the ith shot
    
     % initialize boolean of zeros
    label_RD_ND = boolean(zeros(size(temp_data,1),1)); 
    
    for j = 1:length(parameters)
        
        is_past_threshold = or(temp_data(:,6+j) <= threshold_matrix(j,1),...
            temp_data(:,6+j) >= threshold_matrix(j,2));
        
        label_RD_ND = or(label_RD_ND, is_past_threshold);
        
    end
    
    % if any disruption warning is given, it's a false positive
    if sum(label_RD_ND) > 0 
        n_fp_RD_ND = n_fp_RD_ND + 1;
        list_fp_RD_ND = [list_fp_RD_ND, int32(RD_ND_shots(i))];

    % or else it is a correctly idenitified non-disruption    
    else
        n_correct_RD_ND = n_correct_RD_ND + 1; 
        list_correct_RD_ND = [list_correct_RD_ND, int32(RD_ND_shots(i))];
    end
end

% Now look at the rampdown shots that are "close" to disruption.

n_early_RD_D = 0; % early warnings
n_correct_RD_D = 0; % correctly identified disruptions within t_lower, t_upper
n_late_RD_D = 0; % late warning
n_miss_RD_D = 0; % missed warning. Technically we aren't totally sure
                % because the data doesn't go all the way to 
                % time_until_disrupt == 0
                
% Create lists for shots within the above categories.
list_correct_RD_D = [];
list_early_RD_D = [];
list_miss_RD_D = [];
list_late_RD_D = [];

for i = 1:length(RD_D_shots)
    
    bool_RD_D = ismember(RD_D_data(:,1), RD_D_shots(i));
    
    temp_data = RD_D_data(bool_RD_D, :); % all data for just the ith shot
    
    % initialize boolean of zeros
    label_RD_D = boolean(zeros(size(temp_data,1),1)); 
       
    for j = 1:length(parameters)
        
        is_past_threshold = or(temp_data(:,6+j) <= threshold_matrix(j,1),...
            temp_data(:,6+j) >= threshold_matrix(j,2));
        
        label_RD_D = or(label_RD_D, is_past_threshold);
        
    end
    
    % times before disruption at which a disruption warning is given
    t = temp_data(label_RD_D, 3); 
    
    if isempty(t) % if no disruption warning was given, it was missed!
        n_miss_RD_D = n_miss_RD_D + 1; 
        list_miss_RD_D = [list_miss_RD_D, int32(RD_D_shots(i))];
    else      
        tmax = max(t);
        
        % if the disruption warning is given before t_d s before the disruption
        if tmax > t_upper 
            % we consider it too early (during normal plasma)
            n_early_RD_D = n_early_RD_D +1; 
            list_early_RD_D = [list_early_RD_D, int32(RD_D_shots(i))];
            
        % if the disruption warning is given between t_lower and t_upper
        elseif tmax >= t_lower 
            % it is correctly identified in time
            n_correct_RD_D = n_correct_RD_D + 1;
            list_correct_RD_D = [list_correct_RD_D, int32(RD_D_shots(i))];
         
        % if the warning is within t_lower of the disruption, it's too late    
        else 
            n_late_RD_D = n_late_RD_D + 1;
            list_late_RD_D = [list_late_RD_D, int32(RD_D_shots(i))];
            
        end        
    end
end

disp(['Total flattop disruption shots: ' num2str(length(FD_shots))])
disp(['Correctly identified flattop disruptions: ' num2str(n_correct_FD)])
disp(['Missed disruptions: ' num2str(n_miss_FD)])
disp(['Late warning: ' num2str(n_late_FD)])
disp(['Early warning (false positive): ' num2str(n_early_FD)])
fprintf('\n')

disp(['Total non-disruption shots: ' num2str(length(ND_shots))])
disp(['False positives from ND only: ' num2str(n_fp_ND)])
disp(['Correctly identified non-disruptions: ' num2str(n_correct_ND)])
fprintf('\n')

disp(['Total "non-disruptive" rampdown shots: ' num2str(length(RD_ND_shots))])
disp(['False positives: ' num2str(n_fp_RD_ND)])
disp(['Correctly identified non-disruptive rampdown shots: ' num2str(n_correct_RD_ND)])
fprintf('\n')

disp(['Total "disruptive" rampdown disruption shots: ' num2str(length(RD_D_shots))])
disp(['Correctly identified rampdown disruptions: ' num2str(n_correct_RD_D)])
disp(['Missed rampdown disruptions: ' num2str(n_miss_RD_D)])
disp(['Late warning: ' num2str(n_late_RD_D)])
disp(['Early warning (false positive): ' num2str(n_early_RD_D)])
fprintf('\n')

end