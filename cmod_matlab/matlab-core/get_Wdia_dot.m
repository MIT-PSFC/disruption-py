function [Wdia_dot, Wdia] = get_Wdia_dot(timebase, original_orientation)

% This script gets the diamagnetic energy, Wdia, and calculates its time
% derivative, Wdia_dot, which we are interested in.
%
% Inputs:
%   timebase = array of desired time values
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Outputs:
%   Wdia_dot = time derivative of Wdia, [J/s]
%   Wdia = diamagnetic energy, [J]
%
% Author: Alex Tinguely 2015-09-09

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%mdsopen('mhd', shot);

[wdia, status1] = mdsvalue('\analysis::efit_aeqdsk:wdia'); % actual/measured Wdia
[t, status2] = mdsvalue('dim_of(\analysis::efit_aeqdsk:wdia)'); % time for Wdia

if mod(status1,2)==0 ||mod(status2,2)==0  % if an error occurs
   
    %disp('Error in calling W_dia');
    Wdia_dot = NaN(size(timebase)); % set to NaN
    Wdia = NaN(size(timebase)); % set to NaN
    %return
    
else 

    Wdia = interp1(t, wdia, timebase); % interpolate measured diamagnetic energy over given timebase

    dw = diff(wdia);
    dt = diff(t);
    wdot = dw./dt;

    t2 = t(1:end-1) + 0.5*dt; % create an array at the midpoints of t (one less value than t)

    Wdia_dot = interp1(t2, wdot, timebase);

end

%mdsclose;

% reorient all output arrays to original configuration
Wdia = orient_array(Wdia, original_orientation);
Wdia_dot = orient_array(Wdia_dot, original_orientation);

end