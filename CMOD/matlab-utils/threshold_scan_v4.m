function threshold_scan_v4(shotlist, lower_threshold, upper_threshold, t_upper, t_lower, parameter)

% Alex Tinguely 160411

% The following code was stolen from plot_disruption_data.m.
db = set_database('logbook');

shotlist = unique(sort(shotlist)); % put the shots in order from smallest to largest
                                   % and remove redundant shots if necessary

parameter = lower(strtrim(parameter)); % parameter to be fetched
% just one for this version

% Get data from ALL shots in the shotlist (provided they are in the db)
data = cell2mat(fetch(db, ['select shot, time, time_until_disrupt, ip, ', ... 
    'ip_error, dipprog_dt, ', parameter, ' from disruption_warning ',...
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

% Initialize
n_correct_FD = lower_threshold;
n_early_FD= lower_threshold;
n_late_FD= lower_threshold;
n_miss_FD= lower_threshold;

n_correct_ND= lower_threshold;
n_fp_ND= lower_threshold;

n_correct_RD_ND= lower_threshold;
n_fp_RD_ND= lower_threshold;

n_correct_RD_D= lower_threshold;
n_early_RD_D= lower_threshold;
n_late_RD_D= lower_threshold;
n_miss_RD_D= lower_threshold;

for i = 1:length(lower_threshold)
   
    threshold_matrix = [lower_threshold(i), upper_threshold(i)];
    
    [n_correct_FD(i),n_early_FD(i), n_late_FD(i),n_miss_FD(i), ...
    n_correct_ND(i), n_fp_ND(i), ...
    n_correct_RD_ND(i), n_fp_RD_ND(i), ...
    n_correct_RD_D(i), n_early_RD_D(i),n_late_RD_D(i), n_miss_RD_D(i)] = ...
    threshold_algorithm_for_scan(FD_data, ND_data, RD_data, threshold_matrix, t_lower, t_upper);
    
end

f = figure(6013);
set(f, 'position', [500,500,800,800]);
hold on;
plot(upper_threshold, n_correct_FD, 'k', 'linewidth', 2)
plot(upper_threshold, n_miss_FD + n_late_FD, '--k', 'linewidth', 2)
plot(upper_threshold, n_early_FD, 'r', 'linewidth', 2)
plot(upper_threshold, n_fp_ND, 'b', 'linewidth', 2)
plot(upper_threshold, n_fp_RD_ND + n_early_RD_D + n_late_RD_D + n_miss_RD_D, 'g', 'linewidth', 2)

l =legend('FD: Good Detection', 'FD: Late+Missed', 'FD: Early Warning', ...
    'ND: False Positive', 'RD: FP+EW+LM', ...
    'location', 'best');
set(l, 'fontsize', 14)

ylim([0,length(unique(FD_data(:,1)))])

xlabel(['Threshold value on ' strrep(parameter, '_', ' ')], 'fontsize', 20);

ylabel('Number of shots', 'fontsize', 20);
title([num2str((t_upper)*1e3), 'ms > time until disrupt > ', num2str((t_lower)*1e3), 'ms'], ...
    'fontsize', 20);
set(gca, 'fontsize', 17)

hold off;

end