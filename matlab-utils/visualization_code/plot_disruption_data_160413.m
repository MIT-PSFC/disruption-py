function [d_data, nd_data, d_shots, nd_shots] = plot_disruption_data_160413(shotlist, varargin)

db = set_database('logbook');

shotlist = unique(sort(shotlist)); % put the shots in order from smallest to largest
                                   % and remove redundant shots if necessary

parameter = lower(strtrim(varargin)); % parameters to be fetched

% Get data from ALL shots in the shotlist (provided they are in the db)
data = cell2mat(fetch(db, ['select shot, time, time_until_disrupt, ip, ', ... 
    'ip_error, dipprog_dt, ', parameter{1}, ' from disruption_warning ',...
    'where shot = ', num2str(shotlist(1)),...
    ' order by shot, time']));
% Note the following parameters:
%   data(:,1) = shot
%   data(:,2) = time
%   data(:,3) = time_until_disrupt
%   data(:,4) = ip
%   data(:,5) = ip_error = ip - ip_prog
%   data(:,6) = dipprog_dt

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

data = data(is_flattop,:); % just data during the flattop

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
% non-disruptions don't end before 1s. (This could be upped.)

%nd_bool = zeros(length(nd_shots),1); % initialize 

for i = 1:length(nd_shots)
    
    bool = ismember(nd_data(:,1),nd_shots(i)); % find timeslices for ith nd_shot
    t = nd_data(bool,2); % gets the time array for the ith nd_shot
    
    if max(t) < 1 % if the nd_shot time is less than 1s
        nd_shots(i) = -1; % set the shot to -1 to identify after the loop
    end
end

% Now remove all of the -1's from nd_shots
nd_shots(nd_shots == -1) = [];

% Now we will plot all of our parameters. For disruptions, we will plot as
% a function of time_until_disrupt. For non-disruptions, we will plot as a
% function of time.

% For disruptions

% for each column/parameter of d/nd_data (after the 6 we always call)
for i = 7:size(d_data,2)
   
   % Initialize the number of disruption and non-disruption shots that are
   % actually plotted.   
   n_d = 0;
   n_nd = 0;
    
   % set up the figure
   f = figure(1000 + i - 6);
   set(f, 'position', [500,100,800,1000])
   hold on;
   
   subplot(3,1,1);
   hold on;
   subplot(3,1,2);
   hold on;
   subplot(3,1,3);
   hold on;
   
   for j = 1:length(d_shots)
       
      bool_d = ismember(d_data(:,1), d_shots(j));
      temp_data = d_data(bool_d, :); % all data for just the jth shot
      
      % We should check whether all the data in this shot is just NaNs. If
      % so, then it doesn't contribute to our plot. We can count those that
      % actually are plotted.
      
      % if the number of NaNs is less than the number of data points
      if sum(isnan(temp_data(:,i))) < length(temp_data(:,i))
          n_d = n_d+1;
      end
      
      subplot(3,1,1)
      plot(-temp_data(:,3)*1e3, temp_data(:,i), '-o', 'markersize', 2); 
      xlim([-Inf, 0]);
      
      subplot(3,1,2)
      plot(-temp_data(:,3)*1e3, temp_data(:,i), '-o', 'markersize', 2); 
      xlim([-21, 0]);
       
   end
   
   subplot(3,1,1)
   xlabel('time until disrupt (ms)');
   ylabel([num2str(n_d), ' disruption shots']);
   title(strrep(parameters(i-6), '_', ' '), 'fontsize', 18);
   
   subplot(3,1,2)
   xlabel('time until disrupt (ms)');
   ylabel([num2str(n_d), ' disruption shots']);
 
   for k = 1:length(nd_shots)
       
      bool_nd = ismember(nd_data(:,1), nd_shots(k));
      temp_data2 = nd_data(bool_nd, :); % all data for just the jth shot
      
      % We should check whether all the data in this shot is just NaNs. If
      % so, then it doesn't contribute to our plot. We can count those that
      % actually are plotted.
      
      % if the number of NaNs is less than the number of data points
      if sum(isnan(temp_data2(:,i))) < length(temp_data2(:,i))
          n_nd = n_nd+1;
      end
      
      subplot(3,1,3)
      plot(temp_data2(:,2), temp_data2(:,i), '-o', 'markersize', 2);   
       
   end
   
   subplot(3,1,3)
   xlabel('time (s)');
   ylabel([num2str(n_nd), ' non-disruption shots']);
   hold off;
     
   %saveas(f,['/home/tinguely/Disruptions/Disruption_warning_data_figures/cmod_2015/', parameters{i-6}, '.fig'])
%   saveas(f,['/home/tinguely/Disruptions/Disruption_warning_data_figures/', parameters{i-6}, '.png'])
%   saveas(f,['/home/tinguely/Disruptions/Disruption_warning_data_figures/', parameters{i-6}, '.pdf'])
   %save('/home/tinguely/Disruptions/Disruption_warning_data_figures/cmod_2015/data.mat', 'd_data', 'nd_data', 'parameters')
  
end

end