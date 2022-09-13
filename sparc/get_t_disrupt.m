function time_until_disruption = get_t_disrupt(shot, timebase)

% This function calculates the time until disruption. If the shot isn't a
% disruption, the time until disruption is on the order of 1 million
% seconds which is just unphysical.
%
%  Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   time_until_disruption = t_disrupt - timebase
%
% Author: Alex Tinguely 2015-09-11

[timebase_row, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

shot1 = shot - 0.5; % this range allows us to find ONE shot as specified by shot
shot2 = shot + 0.5;

addpath('/home/granetz/matlab'); % set_database.m is in granetz's matlab folder

db = set_database('logbook');

result = fetch(db, ['select t_disrupt from disruptions where shot between ', ...
    num2str(shot1), ' and ', num2str(shot2)]);
% get the disruption time t_disrupt for the given shot

if isempty(result) % if there is no t_disrupt, i.e. no disruption
    
    t_disrupt = NaN; % sets disruption time to NaN, so time until disruption will also be NaN
    
else

    t_disrupt = result{1};
    
end

time_until_disruption = t_disrupt - timebase_row;

% reorient all output arrays to original configuration
time_until_disruption = orient_array(time_until_disruption, original_orientation);

end