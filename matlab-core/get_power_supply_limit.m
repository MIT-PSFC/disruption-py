function railed_array = get_power_supply_limit(shot, timebase)

% This function is used to determine whether any/all of the power supplies
% have "railed", i.e. the current in the power supply has reached the
% pre-determined limit so it cannot supply any more voltage/power.
% Ideally, if one power supply hits this limit, we want to set the value of
% railed to 1; otherwise, it is zero. So this is a binary system.
%
% One problem that can result is if spacing in timebase is too large, this
% function could totally miss the occurrence of a "railing". An example of
% this is for shot 1120720012: A timebase of 100 elements (~20ms apart)
% catches the railing at ~1.45s, but a timebase of 70 elements, (~30ms
% apart) does NOT.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   railed_array = array of 0 if power supply current is within limits, 
%                  1 if AT LEAST one supply has railed
%
% Author: Alex Tinguely 2015-09-21

[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

% Names of different power supplies to check in the engineering tree.
power_supplies = {'oh1:', 'oh2_u:', 'oh2_l:', 'ef1_u:', 'ef1_l:', 'ef2_u:', ...
                'ef2_l:', 'ef3:u_', 'ef3:l_', 'ef4:l_', 'efc:u_', 'tf:'};

railed_array = zeros(1,length(timebase)); % Initialize
            
mdsopen('engineering', shot);
            
for i = 1:length(power_supplies) % for each power supply
    
   current = mdsvalue(['\engineering::', power_supplies{i}, 'bus_r_cur']); % power supply current [A]
   time = mdsvalue(['dim_of(\engineering::', power_supplies{i}, 'bus_r_cur)']); % time [s] for current
   
   if strcmp(current(1), 'J') || strcmp(time(1), 'J') || length(current) == 1 || length(time) == 1
       
       % if there is an error calling the current or it only has 1 value,
       continue % onto the next power supply, i.e. don't consider this power supply
       % and the railed_array values stay the same
        
   end
   
   current2 = interp1(time, current, timebase); % [A], interpolated current
   
   current_limit = 1e3*mdsvalue(['\engineering::', power_supplies{i}, 'current']); 
   % [A], current limits set by engineering, could be one positive value or
   % array of [negative, postive] limits
   
   for j = 1:length(timebase) % for each time in timebase
      
      % if there are positive and negative limits, and the current at that
      % time exceeds either limit, the power supply has "railed", so we set
      % the value to 1 (otherwise it is zero)
      if length(current_limit) == 2 && (current2(j) < current_limit(1) || current2(j) > current_limit(2))
          railed_array(j) = 1;
      end
      
      % if there is only one positive limit, if the current exceeds this
      % limit or its negative value, set railed_array to 1 at this time
      if length(current_limit) == 1 && (current2(j) < -current_limit(1) || current2(j) > current_limit(1))
          railed_array(j) = 1;
      end
       
   end
   
end

mdsclose;

% reorient all output arrays to original configuration
railed_array = orient_array(railed_array, original_orientation);

end