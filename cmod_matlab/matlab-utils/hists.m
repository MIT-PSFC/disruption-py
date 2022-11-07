% Extract relevant data from the database and define relevant indices

db = set_database('logbook');
get_columns(db,'disruption_warning');
result = fetch(db,['select shot, time, time_until_disrupt, dipprog_dt,'...
    ' p_rad, dprad_dt, li from disruption_warning order by shot,time']);
shot = int64(cell2mat(result(:,1)));
time = cell2mat(result(:,2));
time_until_disrupt = cell2mat(result(:,3));
dipprog_dt = cell2mat(result(:,4));
p_rad = cell2mat(result(:,5));
dprad_dt = cell2mat(result(:,6));
li = cell2mat(result(:,7));
define_indices

% Plot dp_rad/dt vs. time for shots in flattop that didn't disrupt
plot(time(indices_flattop_no_disrupt), ... 
    li(indices_flattop_no_disrupt),'.')
ylim([-4e9,4e9])
xlabel('Time (s)')
ylabel('d_t p_{rad}')
title('All Non-Disruptive Shots in Flattop')

% Plot dp_rad/dt vs. time_until_disrupt for disrupted shots in flattop
figure
plot(-time_until_disrupt(indices_flattop_disrupt_in_flattop), ...
    li(indices_flattop_disrupt_in_flattop),'.')
xlim([-0.03,0])
ylim([-4e9,4e9])
xlabel('Time Before Disrupt (s)')
ylabel('l_{i}')
title('All Disruptive Shots in Flattop')

%{
% Plot same plot but do it shot by shot to see rate of change
figure; hold
shotlist = unique(sort(shot));
for i = 1:length(shotlist)
    shot_indices = find(shot == shotlist(i));
    indices = intersect(shot_indices,indices_flattop_disrupt_in_flattop);
    plot(-time_until_disrupt(indices),p_rad(indices))
end

xlabel('Time Before Disrupt (s)')
ylabel('p_{rad}')
title('All Disruptive Shots in Flattop')
%}    
% Get indices for time periods before disruptions

ind_before20ms = indices_flattop_disrupt_in_flattop(find(...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)>=0.02));

ind_15to20ms = indices_flattop_disrupt_in_flattop(find(...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)>=0.015 & ...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.02));

ind_10to15ms = indices_flattop_disrupt_in_flattop(find(...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)>=0.01 & ...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.015));

ind_20to50ms = indices_flattop_disrupt_in_flattop(find(...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)>=0.05 & ...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.02));

ind_0to20ms = indices_flattop_disrupt_in_flattop(find(...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)>=0 & ...
    time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.020));

bins = [-2e9:1e7:2e9];
x_zoom = [-1e9,1e9];
figure; hold on
histogram(dprad_dt(indices_flattop_no_disrupt), bins, ...
    'normalization','probability');
xlabel('d_t p_{rad}'); xlim(x_zoom);
title('Non-Disruptive Shots in Flattop')

hold off; figure; hold on
histogram(dprad_dt(ind_before20ms),bins, ...
    'FaceColor','r','normalization','probability'); 
histogram(dprad_dt(ind_15to20ms),bins, ...
    'FaceColor','g','normalization','probability');
xlabel('d_t p_{rad}'); xlim(x_zoom);
legend('>20ms before disrupt','15-20ms before','location','northwest')
hold off; figure; hold on
histogram(dprad_dt(ind_10to15ms),bins,'normalization','probability');
histogram(dprad_dt(ind_5to10ms),bins,'normalization','probability');
histogram(dprad_dt(ind_0to5ms),bins,'normalization','probability');
xlabel('d_t p_{rad}'); xlim(x_zoom);
legend('10-15ms before','5-10ms before','0-5ms before', ...
    'location','northwest')
hold off

%{
figure; histogram(dprad_dt(ind_15to20ms),bins); 
title('15-20ms before disruption'); 
xlabel('d_t p_{rad}'); xlim(x_zoom);

figure; histogram(dprad_dt(ind_10to15ms),bins); 
title('10-15ms before disruption');
xlabel('d_t p_{rad}'); xlim(x_zoom);

figure; histogram(dprad_dt(ind_5to10ms),bins); 
title('5-10ms before disruption');
xlabel('d_t p_{rad}'); xlim(x_zoom); 

figure; histogram(dprad_dt(ind_0to5ms),bins); 
title('0-5ms before disruption');
xlabel('d_t p_{rad}'); xlim(x_zoom);

figure; histogram(dprad_dt(indices_flattop_no_disrupt),bins);
title('All Non-Disruptive Shots in Flattop');
xlabel('d_t p_{rad}'); xlim(x_zoom);
%}
n_dis_flattop = size(unique(shot(indices_flattop_disrupt_in_flattop)),1);
n_no_dis_flattop = size(unique(shot(indices_flattop_no_disrupt)),1);