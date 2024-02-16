function [beta_p, beta_p_dot] = get_beta_p(timebase, original_orientation) 

% This script calculates beta_poloidal = beta_p as well as the time
% derivative of beta_p.
%
% Inputs:
%   timebase = array of desired time values
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Outputs:
%   beta_p = beta_poloidal
%   beta_p_dot = time derivative of beta_p
%
% Author: Alex Tinguely 2015-09-09

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%mdsopen('mhd', shot);

betap = mdsvalue('\analysis::efit_aeqdsk:betap'); % actual/measured beta_p
t = mdsvalue('dim_of(\analysis::efit_aeqdsk:betap)'); % time for beta_p

if strcmp(betap(1), 'J') || strcmp(t(1), 'J') || length(betap) == 1 || length(t) == 1
    % if betap,t are not callable or only one point
    beta_p = zeros(1, length(timebase)) + NaN; % return NaN
    beta_p_dot = zeros(1, length(timebase)) + NaN; % return NaN
    return 
end

beta_p = interp1(t, betap, timebase); % interpolate real density over given timebase

dbetap = diff(betap); % difference in betap
dt = diff(t); % difference in time

t2 = t(1:end-1) + 0.5*dt; % create an array at the midpoints of t (one less value than t)

betap_dot = dbetap./dt; % calculate the derivative

beta_p_dot = interp1(t2, betap_dot, timebase); % interpolate over the midpoint array

%mdsclose;

% reorient all output arrays to original configuration
beta_p = orient_array(beta_p, original_orientation);
beta_p_dot = orient_array(beta_p_dot, original_orientation);

end