function [k] = get_elongation(timebase, original_orientation)

% This script calculates elongation of the plasma, k.
%
% Inputs:
%   timebase = array of desired time values
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Outputs:
%   k = elongation
%
% Author: Alex Tinguely 2015-09-09

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%mdsopen('mhd', shot);

K = mdsvalue('\analysis::efit_aeqdsk:eout'); % actual/measured elongation
t = mdsvalue('dim_of(\analysis::efit_aeqdsk:eout)'); % time for elongation

if strcmp(K(1), 'J') || strcmp(t(1), 'J') || length(K) == 1 || length(t) == 1
    % if K,t are not callable or only one point
    k = zeros(1, length(timebase)) + NaN; % return NaN
    return 
end

k = interp1(t, K, timebase); % interpolate elongation over given timebase

% reorient all output arrays to original configuration
k = orient_array(k, original_orientation);

%mdsclose

end