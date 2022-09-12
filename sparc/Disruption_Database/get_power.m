function [p_RAD, p_LH, p_OHM, p_RF, rad_fraction, p_input] = ...
    get_power(shot, timebase, P_ohm, original_orientation)

% NOTE: While this function still works, Bob Granetz has written a better
% version in his directory. That one will be used for the disruption
% warning database, so I will leave this unchanged currently.
% Alex Tinguely, 2015-10-13

% This function gets the input powers -- lower hybrid (LH), ohmic (OHM),
% and radio frequency (RF) -- as well as output radiated power (RAD). If
% any of the requested powers are not callable (there was no RF or LH, for
% instance), then this function sets creates an array of zeros. This should
% really only be a problem for RF and LH, but I've included it for RAD as
% well, just in case. This also calculates the ratio of output radiation to
% input power.
% 
% Inputs:
%   shot = shot number
%   timebase = times to calculate the powers 
%   smooth_span = OPTIONAL argument, the number of elements to
%        average/smooth over in get_P_ohm(). Default is 1, which means no 
%        smoothing. If even, sets to smooth_span-1.
%   original_orientation = 0 if originally column array, 1 if originally row array
%   
% Outputs:
%   p_RAD = radiated power [W]
%   p_LH  = lower hybrid power [W]
%   p_OHM = ohmic power [W]
%   p_RF  = radio frequency power [W]
%   rad_fraction = p_RAD/(p_OHM + p_LH + p_RF), ratio of output radiation to input power
%   p_input = total input power = p_OHM + p_LH + p_RF [W]
%
% Author: Alex Tinguely 15-09-16

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%if nargin == 3 % if smooth_span is provided
%   span = smooth_span;   
%else % if smooth_span is not provided
%   span = 1; % no smoothing at all
%end

% ------------------------------------
% Lower Hybrid
mdsopen('LH', shot);
p_lh = 1e3*mdsvalue('\LH::TOP.RESULTS:NETPOW'); % [W]
t_lh = mdsvalue('dim_of(\LH::TOP.RESULTS:NETPOW)'); % [s]

% if either p_lh or t_lh are not callable, the error *should* start with 
% the word 'Java' and hence letter 'J'  0R the lenght of p_lh,t_lh is 1 (so
% it can't be interpolated)
if strcmp(p_lh(1), 'J') || strcmp(t_lh(1), 'J') || length(p_lh) == 1 || length(t_lh) == 1
   
    p_LH = zeros(1, length(timebase)); % set p_LH to zero
    
else   
    
    p_LH = interp1(t_lh, p_lh, timebase);
    
end

mdsclose;

% ------------------------------------
% RF Power
mdsopen('RF', shot);
p_rf = 1e6*mdsvalue('\rf::rf_power_net'); % [W]
t_rf = mdsvalue('dim_of(\rf::rf_power_net)'); %

if strcmp(p_rf(1), 'J') || strcmp(t_rf(1), 'J') || length(p_rf) == 1 || length(t_rf) == 1
    % if p_rf or t_rf are not callable or only one point
    p_RF = zeros(1, length(timebase));
else   
    p_RF = interp1(t_rf, p_rf, timebase);
end
mdsclose;

% ------------------------------------
% Ohmic power
%p_OHM = get_P_ohm(shot, timebase); %, span); %[W]
p_OHM = P_ohm; % provided to function now

% ------------------------------------
% Radiated power
mdsopen('spectroscopy', shot);
p_rad = mdsvalue('\top.bolometer.twopi_foil'); % [W]
t_rad = mdsvalue('dim_of(\top.bolometer.twopi_foil)'); %

if strcmp(p_rad(1), 'J') || strcmp(t_rad(1), 'J') || length(p_rad) == 1 || length(t_rad) == 1
    % if p_rad or t_rad are not callable or only one point
    p_RAD = zeros(1, length(timebase));
else   
    p_RAD = interp1(t_rad, p_rad, timebase);
end
mdsclose;

% reorient all output arrays to original orientation
p_RAD = orient_array(p_RAD, original_orientation);
p_OHM = orient_array(p_OHM, original_orientation);
p_LH = orient_array(p_LH, original_orientation);
p_RF = orient_array(p_RF, original_orientation);

rad_fraction = p_RAD./(p_OHM + p_LH + p_RF); % = ratio of output to input power

p_input = p_OHM + p_LH + p_RF; % total input power [W]

end