function [n_correct_FD,n_early_FD, n_late_FD,n_miss_FD, ...
          n_correct_ND, n_fp_ND, ...
          n_correct_RD_ND, n_fp_RD_ND, ...
          n_correct_RD_D, n_early_RD_D,n_late_RD_D, n_miss_RD_D] = ...
threshold_algorithm_for_scan(FD_data, ND_data, RD_data, threshold_matrix, t_lower, t_upper)

parameters = 1;

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

end