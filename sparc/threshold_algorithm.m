function [list_d_id, list_d_e, list_d_m, list_d_l, ...
    list_nd_id, list_nd_fp, ...
    list_pd_d, list_pd_nd, list_pd_fp, list_pd_m, list_pd_l] = ...
    threshold_algorithm(shotlist, parameters, threshold_matrix, t_d)

% Here is my attempt at a disruption warning algorithm using thresholds
% I've found while exploring the trends in data as different parameters
% approach a disruption. 
%
% Inputs:
%   shotlist: shots to be tested
%   parameters: string (cell) array of N parameters
%   threshold matrix: 2xN matrix, each row corresponds to one parameter and
%       the 2 columns to the lower and upper thresholds, respectively
%   t_d = time separating pre/disruptions
%
% Alex Tinguely 160407

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
% (should actually be 0, but normal ramp up/down rates are >100kA/s.

ip_prog = data(:,4) - data(:,5); % = ip - ip_error 
is_flattop_1 = abs(data(:,6)) < 6e4; % dipprog_dt < 60kA/s (
is_flattop_2 = abs(ip_prog) > 1e5; % ip_prog > 100kA
is_flattop = and(is_flattop_1, is_flattop_2);

data = data(is_flattop,:); % just data during the flattop for both non/disruptions

% Now we want to split up our data into disruptions and non-disruptions.
% First, we check if the time_until_disrupt is NaN, which is indicative of
% a non-disruption. 

t_is_nan = isnan(data(:,3)); % time_until_disrupt is not NaN

nd_data = data(t_is_nan,:);
d_data = data(~t_is_nan,:);

% We only want disruptions that actually occur during the flattop.
% Currently we have flattop data, but there is the possibility that
% time_until_disrupt = 0 occurs after the flattop. We want to get rid of
% these shots.

is_t_until_disrupt_0 = d_data(:,3) == 0; % boolean if time_until_disrupt == 0
% find shot numbers where disrupton occurs during flattop
shots_t0 = unique(d_data(is_t_until_disrupt_0, 1));
% get a boolean array of the timeslices for these shots
is_shot_with_t0 = ismember(d_data(:,1),shots_t0);
D_data = d_data(is_shot_with_t0,:);

PD_data = d_data(~is_shot_with_t0,:);
ND_data = nd_data;

% Get the shots in each
ND_shots = sort(unique(ND_data(:,1)));
PD_shots = sort(unique(PD_data(:,1)));
D_shots = sort(unique(D_data(:,1)));

% Now we will walk through each shot and count the number of false
% positives, early warnings, etc.

% initialize the number of...
n_early = 0; % early warnings
n_d_id = 0; % correctly identified disruptions
n_late = 0; % late warning
n_miss = 0; % missed disruption
n_fp = 0; % false positive
n_nd_id = 0; % correctly identified non-disruptions

list_d_id = [];
list_d_e = [];
list_d_m = [];
list_d_l = [];

% First consider the disruption shots.
for i = 1:length(D_shots)
    
    % get boolean array for all timeslices of the ith disruption test shot.
    % This uses d_test_data from above which inludes all timeslices for the
    % disruption (even t_until_disrupt > 20ms)
    bool_D = ismember(D_data(:,1), D_shots(i));
    
    temp_data = D_data(bool_D, :); % all data for just the ith shot
    
    % initialize boolean of zeros
    label_D = boolean(zeros(size(temp_data,1),1)); 
    
    % For each of the parameters (starting at column 7), we see if the data
    % is past the thresholds as provided in threshold_matrix. Then we use
    % the or() function to compare with the other parameters. This
    % "label_D" then indicates whether the datapoint was considered to be a
    % disruption (1 = yes, 0 = no).
    
    for j = 1:length(parameters)
        
        is_past_threshold = or(temp_data(:,6+j)<=threshold_matrix(j,1),...
            temp_data(:,6+j)>=threshold_matrix(j,2));
        
        label_D = or(label_D, is_past_threshold);
        
    end
    
    t = temp_data(label_D, 3); % times before disruption at which a disruption warning is given
    
    if isempty(t) % if no disruption warning was given, it was missed!
        n_miss = n_miss + 1; 
        list_d_m = [list_d_m, int32(D_shots(i))];
    else
        %figure(800+i); plot(-temp_data(:,3)*1e3, label_D); % plot -time_until_disrupt vs. label
        %xlabel('time until disrupt (ms)');
        tmax = max(t);
        if tmax > t_d % if the disruption warning is given before t_d s before the disruption
            n_early = n_early +1; % we consider it too early (during normal plasma)
            
            list_d_e = [list_d_e, int32(D_shots(i))];
            %title(['Early warning: ', num2str(d_test_shots(i))]);
            %xlim([-tmax*1e3 - 5,0]);
        elseif tmax >= 0.0095 % if the disruption is given between t_d and 9.5ms,
            n_d_id = n_d_id + 1; % it is correctly identified in time
            list_d_id = [list_d_id, int32(D_shots(i))];
            %title(['Correctly identified disruption: ', num2str(d_test_shots(i))]);
            %xlim([-25,0]);
        else % if the warning is within 10ms of the disruption, it's too late
            n_late = n_late + 1;
            list_d_l = [list_d_l, int32(D_shots(i))];
            %title(['Late warning: ', num2str(d_test_shots(i))]);
            %xlim([-25,0]);
        end        
    end
end

list_nd_id = [];
list_nd_fp = [];


% Next look at the non-disruption shots
for i = 1:length(ND_shots)
    
    % get indices for all timeslices of the ith non-disruption test shot
    %bool_NDtest = ismember(nd_data(:,1), nd_test_shots(i));
    
    % This uses ND_data from above, but still have to select out the
    % data for each test shot.
    bool_ND2 = ismember(ND_data(:,1), ND_shots(i));
    temp_data = ND_data(bool_ND2, :); % all data for just the ith shot
    
     % initialize boolean of zeros
    label_ND = boolean(zeros(size(temp_data,1),1)); 
    
    for j = 1:length(parameters)
        
        is_past_threshold = or(temp_data(:,6+j)<=threshold_matrix(j,1),...
            temp_data(:,6+j)>=threshold_matrix(j,2));
        
        label_ND = or(label_ND, is_past_threshold);
        
    end
    
    if sum(label_ND) > 0 % if any disruption warning was given, it's a false positive
        n_fp = n_fp + 1;
        list_nd_fp = [list_nd_fp, int32(ND_shots(i))];
        %figure(700+i); plot(temp_data(:,2), label_ND); % plot time vs. label
        %title(['False positive: ', num2str(nd_test_shots(i))]);
        %xlabel('time (s)');
    else
        n_nd_id = n_nd_id + 1; % correctly idenitified non-disruptions
        list_nd_id = [list_nd_id, int32(ND_shots(i))];
    end
end

% Now we look at shots that are like "pre" disruptions, where they disrupt
% *after* the flattop. Thus, for the most part, they are just like
% non-disruptions. However, if the disruption occurs relatively soon after
% the flattop end (within t_d s), then we can sort of treat this like a
% disruption. If this is the case, can we even predict these disruptions?

n_pd_nd = 0; % correctly idenitfy pre-disruption as non-disruptive
n_miss_pd = 0;
n_fp_pd = 0; % both false positive and early
n_pd_d = 0; % correctly idenitfy pre-disruption as disruptive
n_late_pd = 0;

list_pd_d = [];
list_pd_nd = [];
list_pd_fp = [];
list_pd_m = [];
list_pd_l = [];

% Finally consider the "pre"-disruption shots.
for i = 1:length(PD_shots)
    
    % get boolean array for all timeslices of the ith disruption test shot.
    % This uses PD_data from above which inludes all timeslices for the
    % disruption (even t_until_disrupt > 20ms)
    bool_PD = ismember(PD_data(:,1), PD_shots(i));
    
    temp_data = PD_data(bool_PD, :); % all data for just the ith shot
    
    % initialize boolean of zeros
    label_PD = boolean(zeros(size(temp_data,1),1)); 
    
    % For each of the parameters (starting at column 7), we see if the data
    % is past the thresholds as provided in threshold_matrix. Then we use
    % the or() function to compare with the other parameters. This
    % "label_PD" then indicates whether the datapoint was considered to be a
    % disruption (1 = yes, 0 = no).
    
    for j = 1:length(parameters)
        
        is_past_threshold = or(temp_data(:,6+j)<=threshold_matrix(j,1),...
            temp_data(:,6+j)>=threshold_matrix(j,2));
        
        label_PD = or(label_PD, is_past_threshold);
        
    end
    
    t = temp_data(label_PD, 3); % times before disruption at which a disruption warning is given
    
    % We are also interested in whether there occur any times in this
    % "pre"-disruption that are less than the given t_d s. If there are,
    % then this means that this pre-disruption shot disrupted within t_d s
    % of the end of the flattop.
    t_leq_t_d = temp_data(:,3) <= t_d; % leq = less than or equal to
    
    % if no disruption warning was given within t_d s, it was missed!
    if isempty(t) && sum(t_leq_t_d) > 0 
        n_miss_pd = n_miss_pd + 1;
        list_pd_m = [list_pd_m, int32(PD_shots(i))];
        
    % if no warning was given, but this pd is just like a nd
    elseif isempty(t) && sum(t_leq_t_d) == 0
        n_pd_nd = n_pd_nd + 1;
        list_pd_nd = [list_pd_nd, int32(PD_shots(i))];
    else
        %figure(800+i); plot(-temp_data(:,3)*1e3, label_D); % plot -time_until_disrupt vs. label
        %xlabel('time until disrupt (ms)');
        tmax = max(t);
        if tmax > t_d % if the disruption warning is given before t_d s before the disruption
            n_fp_pd = n_fp_pd +1; % we consider it a false-positive/too early
            list_pd_fp = [list_pd_fp, int32(PD_shots(i))];
            %title(['Early warning: ', num2str(d_test_shots(i))]);
            %xlim([-tmax*1e3 - 5,0]);
        elseif tmax >= 0.0095 % if the disruption is given before 9.5ms,
            n_pd_d = n_pd_d + 1; % it is correctly identified
            
            list_pd_d = [list_pd_d, int32(PD_shots(i))];
            
            %title(['Correctly identified disruption: ', num2str(d_test_shots(i))]);
            %xlim([-25,0]);
        else % if the warning is within 10ms of the disruption, it's too late
            n_late_pd = n_late_pd + 1;
            list_pd_l = [list_pd_l, int32(PD_shots(i))];
            %title(['Late warning: ', num2str(d_test_shots(i))]);
            %xlim([-25,0]);
        end        
    end
end


disp(['Correctly identified disruptions (10-20ms): ' num2str(n_d_id)])
disp(['Missed disruptions: ' num2str(n_miss)])
disp(['Late warning: ' num2str(n_late)])
disp(['Early warning (false positive): ' num2str(n_early)])
disp(['Total disruption shots: ' num2str(length(D_shots))])
fprintf('\n')
disp(['False positives from ND only: ' num2str(n_fp)])
disp(['Correctly identified non-disruptions: ' num2str(n_nd_id)])
disp(['Total non-disruption shots: ' num2str(length(ND_shots))])
fprintf('\n')
disp(['Correctly identify pre-disruptions as disruptions (10-20ms): ' num2str(n_pd_d)])
disp(['Correctly identify pre-disruptions as non-disruptions: ' num2str(n_pd_nd)])
disp(['"Missed" pre-disruptions: ' num2str(n_miss_pd)])
disp(['"Late" warning: ' num2str(n_late_pd)])
disp(['False positive for pre-disruption: ' num2str(n_fp_pd)])
disp(['Total pre-disruption shots: ' num2str(length(PD_shots))])

end