function [Li, Li_dot] = get_Li(timebase, original_orientation)

% This script gets the internal inductance, Li, and calculates its time
% derivative, Li_dot.
%
% Inputs:
%   timebase = array of desired time values
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Outputs:
%   Li = internal inductance
%   Li_dot = time derivative of Li
%
% Author: Alex Tinguely 2015-09-09

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%mdsopen('mhd', shot);

li = mdsvalue('\analysis::efit_aeqdsk:li'); % actual/measured Li
t = mdsvalue('dim_of(\analysis::efit_aeqdsk:li)'); % time for Li

if strcmp(li(1), 'J') || strcmp(t(1), 'J') || length(li) == 1 || length(t) == 1
    % if li,t are not callable or only one point
    Li = zeros(1, length(timebase)) + NaN; % return NaN
    Li_dot = zeros(1, length(timebase)) + NaN; % return NaN
    return 
end

Li = interp1(t, li, timebase); % interpolate measured inductance over given timebase

dli = diff(li);
dt = diff(t);

t2 = t(1:end-1) + 0.5*dt; % create an array at the midpoints of t (one less value than t)

li_dot = dli./dt;

Li_dot = interp1(t2, li_dot, timebase);

%mdsclose;

% reorient all output arrays to original configuration
Li = orient_array(Li, original_orientation);
Li_dot = orient_array(Li_dot, original_orientation);

end