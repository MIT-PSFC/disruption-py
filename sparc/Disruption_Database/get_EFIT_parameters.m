function [beta_N, beta_p, beta_p_dot, k, Upper, Lower, Li, Li_dot, ...
    p_RAD, p_LH, p_OHM, p_RF, rad_fraction, p_input, q0, qstar, q95, ...
    V_loop, Wdia_dot, Wdia] = get_EFIT_parameters(shot, timebase)

% This function calculates *most* of the EFIT parameters as well as a few
% non-EFIT parameters. The reasoning behind this function is that we need
% to pull out many EFIT parameters, but it wastes time for each individual
% function to mdsopen 'mhd' and mdsclose it. Therefore, we mdsopen once at
% the beginning, calculate all the values, and mdsclose it.
%
% The only non-EFIT parameters are those in get_power, which get data from
% LH tree, etc. The reason that I include it here is that P_ohm is
% calculated from EFIT and is required in the p_input and rad_fraction
% calculations of get_power(). Because this is at the end, the
% mdsopening/mdsclosing of those trees doesn't affect the previous EFIT
% data.
%
% Inputs:
%   shot = shot number
%   timebase = array of times at which to evaluate the parameters
%
% Outputs:
%   beta_N = normalized beta
%   beta_p = beta_poloidal
%   beta_p_dot = time derivative of beta_p [1/s]
%   k = elongation
%   Upper = upper gap [m]
%   Lower = lower gap [m]
%   Li = internal inductance [H]
%   Li_dot = time derivative of Li [H/s]
%   p_OHM = ohmic power [W]
%   q0 = safety factor at center
%   qstar = q* = cylindrical safety factor
%   q95 = edge safety factor
%   V_loop = loop voltage [V]
%   Wdia_dot = time derivative of Wdia, [J/s]
%   Wdia = diamagnetic energy, [J]
%   p_RAD = radiated power [W]
%   p_LH  = lower hybrid power [W]
%   p_RF  = radio frequency power [W]
%   rad_fraction = p_RAD/(p_OHM + p_LH + p_RF), ratio of output radiation to input power
%   p_input = total input power = p_OHM + p_LH + p_RF [W]
%
% Author: Alex Tinguely 2015-10-07 (Dan Brunner's 30th Birthday!)

% make the timebase a row array, easier to do calculations this way
[timebase, original_orientation] = make_row_array(timebase);
% original_orientation is 0 if originally a column array, 1 if originally row array

mdsopen('mhd', shot); % to get EFIT parameters

beta_N = get_beta_N(timebase, original_orientation);
[beta_p, beta_p_dot] = get_beta_p(timebase, original_orientation);
k = get_elongation(timebase, original_orientation);
[Upper, Lower] = get_gaps(timebase, original_orientation);
[Li, Li_dot] = get_Li(timebase, original_orientation);
p_OHM = get_P_ohm(timebase, original_orientation);
[q0, qstar, q95] = get_safety_factor(timebase, original_orientation);
V_loop = get_V_loop(timebase, original_orientation);
[Wdia_dot, Wdia] = get_Wdia_dot(timebase, original_orientation);

mdsclose; % get out of EFIT

[p_RAD, p_LH, ~, p_RF, rad_fraction, p_input] = get_power(shot, timebase, p_OHM, original_orientation);

end