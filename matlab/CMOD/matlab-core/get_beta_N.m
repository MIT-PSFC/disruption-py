function [beta_N] = get_beta_N(timebase, original_orientation)

% This script calculates beta_N = normalized beta
%
% Inputs:
%   timebase = array of desired time values
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Outputs:
%   beta_N = beta_Normalized
%
% Author: Alex Tinguely 2015-09-09

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%mdsopen('mhd', shot);

betaN = mdsvalue('\analysis::efit_aeqdsk:betan'); % actual/measured beta_N
t = mdsvalue('dim_of(\analysis::efit_aeqdsk:betan)'); % time for beta_N

if strcmp(betaN(1), 'J') || strcmp(t(1), 'J') || length(betaN) == 1 || length(t) == 1
    % if betaN,t are not callable or only one point
    beta_N = zeros(1, length(timebase)) + NaN; % return NaN
    return 
end

beta_N = interp1(t, betaN, timebase); % interpolate beta_N over given timebase

%mdsclose;

% reorient all output arrays to original configuration
beta_N = orient_array(beta_N, original_orientation);

end