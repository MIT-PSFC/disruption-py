function [D_data, PD_data, ND_data] = ...
    plot_disruption_histogram(shotlist, bin_width, xlimits, edges, t_d, parameter)

% This script plots a histogram of our disruption warning database data at
% each timeslice from 20ms to 0ms before a disruption. Because some EFIT
% times can lie between each 1ms interval, the plots are for intervals t>=
% time > t-1 ms. These are then recorded in a movie to see trends in time.
% The histograms for non/disruptions are overlayed, both normalized so all
% points add to one. Note that all times before 20ms before a disruption
% are considered non-disruptions.

% Inputs:
%   shotlist = list of shots to get data
%   bin_width = width of bins in the histogram
%   xlimits = array limits for the plot in x
%   edges = array of edges between which to calculate the number of counts
%   t = time (s) before disruption used to separate pre/disruptive data
%   parameter = string, database parameter to be tested (currently can only
%   do one)

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

%return

% Now we will plot histograms of our data for each parameter of interest.
% First we will plot a histogram for all non/pre-disruptive data (over all
% times available) so this doesn't change. Then, for each ms before the
% disruption, we will update the histogram for disruptive data overlayed on
% the non-disruptive data.

v = VideoWriter(['/home/tinguely/Disruptions/Disruption_warning_data_figures/', ...
    'cmod_2015/histogram_', parameter,'.avi']);

open(v);

h_cell = cell(20,1); % initialize

for i = 1:20 % ms
%for j = 11; % this will just be the case of 10ms >=t > 9ms
    
    % open a (temporary) figure
    f = figure(2000+i);

    % The histogram for ND_data is over all times
    h_nd = histogram(ND_data(:,7), 'facecolor', 'b', 'facealpha', 0.5);
    h_nd.Normalization = 'probability'; % normalize to probability
    h_nd.BinWidth = bin_width; % set bin width
    
    hold on;
    
    % Now get the number of counts between the specified "edges" and
    % compute the probability.
    
    [N_nd] = histcounts(ND_data(:,7), edges);
    P_nd = N_nd./sum(N_nd);
    
    % The histogram for PD_data (pre-disruptive) is over all times
    h_pd = histogram(PD_data(:,7), 'facecolor', 'r', 'facealpha', 0.5);
    h_pd.Normalization = 'probability'; % normalize to probability
    h_pd.BinWidth = bin_width; % set bin width
    
    % Now get the number of counts between the specified "edges" and
    % compute the probability.
    
    [N_pd] = histcounts(PD_data(:,7), edges);
    P_pd = N_pd./sum(N_pd);
        
    % We also want to know the percentage of non-disruption *shots* that are
    % false positives or correctly identified. False positive means that
    % the data points are outside our edge "limits" so they are beyond the
    % threshold.
    
    is_outside_edges = or(ND_data(:,7)<= edges(2),ND_data(:,7)>= edges(3));
    false_pos = length(unique(ND_data(is_outside_edges,1)))/length(unique(ND_data(:,1)));
    
    % Find the percentage of pre-disruption shots that qualify as "early".
    % Again, this is data outside the edges, but before t_d ms before a
    % disruption, so it would be an early warning.
    
    is_outside_edges = or(PD_data(:,7)<= edges(2),PD_data(:,7)>= edges(3));
    early = length(unique(PD_data(is_outside_edges,1)))/length(unique(PD_data(:,1)));

    % Find D_data where time_until_disrupt is between 20-i and 21-i ms,
    % i.e. we are stepping through each timeslice from 20ms to 0ms
    % before a disruption

    bool3 = and(D_data(:,3) < (21.5-i)*1e-3, D_data(:,3) > (20.5-i)*1e-3);
    D_temp = D_data(bool3,:); % make a temp D_data matrix for this timeslice
    
    h_d = histogram(D_temp(:,7), 'facecolor', 'g', 'facealpha', 0.7);
    h_d.Normalization = 'probability'; % normalize to probability
    h_d.BinWidth = bin_width;

    h_cell{i} = h_d;

    % Now get the number of counts between the specified "edges" and
    % compute the probability.
    
    [N_d] = histcounts(D_temp(:,7), edges);
    P_d = N_d./sum(N_d);
    
    % Find the percentage of disruptions correctly identified and missed.
    
    is_outside_edges = or(D_temp(:,7)<= edges(2),D_temp(:,7)>= edges(3));
    correct = length(unique(D_temp(is_outside_edges,1)))/length(unique(D_temp(:,1)));
    missed = 1-correct;
    
    % Now we plot the edges. We expect an array of length 4 given in the
    % form [-Inf, a, b, Inf], so we only want to plot vertical lines at a
    % and b.
    
    ax = gca;
    ylimits = ax.YLim;
    plot([edges(2), edges(2)], [ylimits(1), ylimits(2)], 'r', 'LineWidth', 2);
    plot([edges(3), edges(3)], [ylimits(1), ylimits(2)], 'r', 'LineWidth', 2);
    
    % Add a legend.
    
%     legend([num2str(length(ND_data(:,i))), ' ND timeslices'],...
%         [num2str(length(PD_data(:,i))), ' PD timeslices'], ...
%         [num2str(length(D_temp(:,i))), ' D timeslices'], ...
%         num2str(edges(2)), num2str(edges(3)), 'Location', 'southoutside');
    
    legend([num2str(length(unique(ND_data(:,1)))), ' Non-Disruption shots'],...
        [num2str(length(unique(PD_data(:,1)))), ' All Disruption shots'], ...
        [num2str(length(unique(D_temp(:,1)))), ' Flattop disruption shots'], ...
        num2str(edges(2)), num2str(edges(3)), 'Location', 'southoutside');


    xlabel(strrep(parameter, '_', ' '), 'fontsize', 18);
    ylabel('probability');
    title([num2str(21.5-i), 'ms > time until disrupt > ', num2str(20.5-i), 'ms'], ...
        'fontsize', 18);

    xlim(xlimits);
    ylim(ylimits);

    % Add text to plot with the percentage of counts inside and outside the
    % edge limits.
    
%     text(edges(2), ylimits(2)*0.95, num2str(P_nd(2)), 'Color', 'b');
%     text(edges(2), ylimits(2)*0.9, num2str(P_d(2)), 'Color', 'r');
%     text(edges(3), ylimits(2)*0.85, num2str(P_nd(1)+P_nd(3)), 'Color', 'b');
%     text(edges(3), ylimits(2)*0.8, num2str(P_d(1)+P_d(3)), 'Color', 'r');
    
    text(ax.XLim(1), ylimits(2)*0.95, ['false positive: ' num2str(false_pos)], 'Color', 'b');
    text(ax.XLim(1), ylimits(2)*0.9, ['early warning: ' num2str(early)], 'Color', 'r');
    text(ax.XLim(1), ylimits(2)*0.85, ['missed warning: ' num2str(missed)], 'Color', 'g');
    text(ax.XLim(1), ylimits(2)*0.8, ['correct id disruption: ' num2str(correct)], 'Color', 'k');    
    text(ax.XLim(1), ylimits(2)*0.75, ['disruption data within ' num2str(t_d) 'ms'], 'Color', 'c');  
    
    set(f, 'position', [500,500,700,700]);
    frame = getframe(f);        
    writeVideo(v,frame);

%     display(['fp: ' num2str(false_pos)])
%     display(['ew: ' num2str(early)])
%     display(['cd: ' num2str(correct)])
%     display(['m: ' num2str(missed)])
        
end

%close all the figures. For some reason the code only works this way.
for i = 1:20
    close(figure(2000+i));
end

close all
close(v);

implay(['/home/tinguely/Disruptions/Disruption_warning_data_figures/', ...
    'cmod_2015/histogram_', parameter,'.avi']);
    
%end

end