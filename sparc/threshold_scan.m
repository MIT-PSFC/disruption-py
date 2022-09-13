function threshold_scan(shotlist, lower_threshold, upper_threshold, t_d, T, parameter)

% Alex Tinguely 160406

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

% We want to know all of our non/disruption shots from d/nd_data

d_shots = sort(unique(d_data(:,1))); % all unique disruption shots
nd_shots = sort(unique(nd_data(:,1))); % all unique non-disruption shots

% We also want to be careful of "non-disruptions" that last for less
% (approximately) the full two seconds. We have now restricted our data to
% those timeslices during the flattop portion (~1-1.5s). Let's make sure our
% *non-disruptions* don't end before 1s. (This could be upped.)

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

% For our purposes, we only consider datapoints within t s of a disruption
% to actually be disruptive. Therefore, all disruptive timeslices more than
% t s away could be called "pre-disruptive" and we could add them to our
% non-disruption data. Note that usually, t = 20ms, but we can specify that
% when we call the function.

is_less_than_t = d_data(:,3) <= t_d; % t_until_disrupt <= t (s)
D_data = d_data(is_less_than_t, :);
PD_data = d_data(~is_less_than_t, :);
ND_data = nd_data;
%ND_data = [nd_data; PD_data];

% We only want disruptions that actually occur during the flattop.
% Currently we have flattop data, but there is the possibility that
% time_until_disrupt = 0 occurs after the flattop. We want to get rid of
% these shots.

is_t_until_disrupt_0 = D_data(:,3) == 0; % boolean if time_until_disrupt == 0
% find shot numbers where disrupton occurs during flattop
shots_t0 = unique(D_data(is_t_until_disrupt_0, 1));
% get a boolean array of the timeslices for these shots
is_shot_with_t0 = ismember(D_data(:,1),shots_t0);
D_data = D_data(is_shot_with_t0,:);

% Also get rid of negative time_until_disrupt data.
D_data(D_data(:,3) < 0,:) = [];

% Look at the timeslice around 10ms before a disruption only.

% We also want to know the percentage of non-disruption *shots* that are
% false positives or correctly identified. False positive means that
% the data points are outside our edge "limits" so they are beyond the
% threshold.

% Initialize
false_pos = lower_threshold;
early = lower_threshold;
correct = lower_threshold;
missed = lower_threshold;

for i = 1:length(lower_threshold)

    is_past_threshold = or(ND_data(:,7)<= lower_threshold(i),ND_data(:,7)>= upper_threshold(i));
    %false_pos(i) = length(unique(ND_data(is_past_threshold,1)))/length(unique(ND_data(:,1)));
    false_pos(i) = length(unique(ND_data(is_past_threshold,1)));
    
    % Find the percentage of pre-disruption shots that qualify as "early".
    % Again, this is data outside the edges, but before t_d ms before a
    % disruption, so it would be an early warning.

    is_past_threshold = or(PD_data(:,7)<= lower_threshold(i),PD_data(:,7)>= upper_threshold(i));
    %early(i) = length(unique(PD_data(is_past_threshold,1)))/length(unique(PD_data(:,1)));
    early(i) = length(unique(PD_data(is_past_threshold,1)));
    
    % Find D_data where time_until_disrupt is between t +- 0.005 s

    bool3 = and(D_data(:,3) < (T+0.0005), D_data(:,3) > (T - 0.0005));
    D_temp = D_data(bool3,:); % make a temp D_data matrix for this timeslice

    is_past_threshold = or(D_temp(:,7)<= lower_threshold(i),D_temp(:,7)>= upper_threshold(i));
    %correct(i) = length(unique(D_temp(is_past_threshold,1)))/length(unique(D_temp(:,1)));   
    correct(i) = length(unique(D_temp(is_past_threshold,1)));   
    missed(i) = length(unique(D_temp(:,1))) - length(unique(D_temp(is_past_threshold,1))); 
    
end

%missed = 1-correct;


f = figure(3000);
set(f, 'position', [500,500,700,700]);
hold on;
plot(upper_threshold, correct, 'k', 'linewidth', 2)
plot(upper_threshold, missed, '--k')
plot(upper_threshold, false_pos, 'b', 'linewidth', 2)
plot(upper_threshold, early, 'r', 'linewidth', 2)

ylim([0,correct(1)])

legend('Correctly ID Disruption', 'Missed/Late Disruption', 'False Positive', 'Early Warning', ...
    'location', 'southoutside')
legend('Correctly ID Disruption', 'Missed/Late Disruption', 'False Positive', 'Early Warning', ...
    'location', 'best')

xlabel(['Threshold Value on ' strrep(parameter, '_', ' ')], 'fontsize', 18);
%ylabel('Percentage');
ylabel('Number of shots');
title([num2str((T + 0.0005)*1e3), 'ms > time until disrupt > ', num2str((T - 0.0005)*1e3), 'ms'], ...
    'fontsize', 18);

ax = gca;
% text(ax.XLim(2)*.55, ax.YLim(2)*0.95, ['Non-Disruption shots: ' num2str(length(unique(ND_data(:,1))))], 'Color', 'b');
% text(ax.XLim(2)*.55, ax.YLim(2)*0.9, ['All Disruption shots: ' num2str(length(unique(PD_data(:,1))))], 'Color', 'r');
% text(ax.XLim(2)*.55, ax.YLim(2)*0.85, ['Flattop Disruption shots: ' num2str(length(unique(D_data(:,1))))], 'Color', 'k');
% text(ax.XLim(2)*.55, ax.YLim(2)*0.8, ['Disruption data within ' num2str(t_d*1e3) 'ms'], 'Color', 'g');

hold off;

end
