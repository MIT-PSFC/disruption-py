function [Upper, Lower] = get_gaps(timebase, original_orientation)

% This script calculates upper and lower gaps between the plasma and the
% upper and lower walls (respectively).
%
% Inputs:
%   timebase = array of desired time values
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Outputs:
%   Upper = upper gap (m)
%   Lower = lower gap (m)
%
% Author: Alex Tinguely 2015-09-09

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%mdsopen('mhd', shot);

upper = mdsvalue('\analysis::efit_aeqdsk:otop')*1e-2; % [m] actual/measured upper gap
tu = mdsvalue('dim_of(\analysis::efit_aeqdsk:otop)'); % time

lower = mdsvalue('\analysis::efit_aeqdsk:obott')*1e-2; % [m] actual/measured lower gap
tl = mdsvalue('dim_of(\analysis::efit_aeqdsk:obott)'); % time

if strcmp(upper(1), 'J') || strcmp(tu(1), 'J') || length(upper) == 1 || length(tu) == 1
    % if upper,ut are not callable or only one point
    Upper = zeros(1, length(timebase)) + NaN; % return NaN
else
    Upper = interp1(tu, upper, timebase); % interpolate upper gap over given timebase 
end

if strcmp(lower(1), 'J') || strcmp(tl(1), 'J') || length(lower) == 1 || length(tl) == 1
    % if lower,tl are not callable or only one point
    Lower = zeros(1, length(timebase)) + NaN; % return NaN 
else    
    Lower = interp1(tl, lower, timebase); % interpolate lower gap over given timebase
end

% reorient all output arrays to original configuration
Upper = orient_array(Upper, original_orientation);
Lower = orient_array(Lower, original_orientation);

%mdsclose

end