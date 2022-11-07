function [V_loop] = get_V_loop(timebase, original_orientation)

% This script calculates loop voltage.
%
% Inputs:
%   timebase = array of desired time values
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Outputs:
%   V_loop = loop voltage
%
% Author: Alex Tinguely 2015-09-09

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%mdsopen('mhd', shot);

[vloop, status1] = mdsvalue('\analysis::top.mflux:v0'); % actual/measured 
[t, status2] = mdsvalue('dim_of(\analysis::top.mflux:v0)'); % time

if mod(status1,2)==0 || mod(status2,2)==0 % if an error occurs
   
    %disp('Error in calling V_loop');
    V_loop = NaN(size(timebase)); % set to NaN
    %return
    
else

    V_loop = interp1(t, vloop, timebase); % interpolate v_loop over given timebase

end
%mdsclose;

% reorient all output arrays to original configuration
V_loop = orient_array(V_loop, original_orientation);

end