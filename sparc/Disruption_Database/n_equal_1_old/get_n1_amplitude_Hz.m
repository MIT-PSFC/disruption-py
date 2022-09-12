function [amp_array] = get_n1_amplitude_Hz(shot, timebase, num_times)

% This function finds the best fit of the n = 1  MHD mode around the torus for
% five different toroidal locations (Mirnov coils) using fluctuations in
% the TIME DERIVATIVE of the poloidal magnetic field (T/s). This is based off of
% work done by Geoff Olynyk in his thesis. Then this is normalized by the
% toroidal magnetic field, so the resulting n=1 amplitude ratio is in Hz.
%
% The idea is that there are constantly fluctuations in the poloidal
% magnetic field over time, and these fluctuations vary in different
% locations around the tokamak. If at each point in time, we can look at
% the fluctuations at different toroidal locations, we can determine the
% SPATIAL variation/fluctuation in the magnetic field. This is the n MHD
% mode number. We want to look for an increase in the n = 1 mode spatial
% fluctuation amplitude preceding a disruption, which will help warn us. 
%
% NOTE: This code was adapted from fit_n_mode.m in
% /home/tinguely/Disruptions/Code. Also, this takes about 30-75 s to
% compute depending on the input timebase and num_times.
%
% Author: Alex Tinguely 2015-09-14
%
% Inputs:
%   shot = shot number
%   timebase = array of requested times
%   num_times = OPTIONAL ARGUMENT, number of time slices to evaluate the n=1 mode
%               fit. A greater number of divisions (large num_times) means
%               a better fit. If no argument is provided, it will default
%               to the number of elements in timebase.
%
% Outputs:
%   amp_array = array of n=1 mode amplitudes in dB_poloidal/dt = (dB_p,n=1)/B_T
%  

if nargin == 3 % if num_times is provided
   
    time = linspace(timebase(1), timebase(end), num_times);
    % make an array of times between the start and stop times of timebase
    % but with num_times number of elements
    
else % if num_times is not provided, just use timebase
    
    time = timebase; 

end

[time, original_orientation] = make_row_array(time); % makes timebase into row array
% easier to do calculations this way

format long % increase the precision of the numbers in the mdsvalues

bp_str = {'_AB_TOP', '_BC_TOP', '_EF_TOP', '_KA_TOP', '06_GHK'}; % 5 Mirnov coil names
phi = [349.8, 300.1, 190.4, 15.2, 128.1]; % 5 toroidal locations (degrees)
rad = phi*pi/180; % in radians

lb = [0, -pi]; % lower bound for fit [amplitude, phase], don't want negative amplitudes
ub = [Inf, pi]; % upper bound for fit [amplitude, phase]
option = optimset('TolFun', 0, 'MaxFunEvals', 1e8, 'Display', 'off');
% in the options, increase the number of function evaluations and suppress
% the display of text after every iteration
    
mdsopen('magnetics', shot);

bp = cell(1,5); % initialize cell array for all 5 locations
b_matrix = zeros(length(time), 5);

for j = 1:5 % for each toroidal location

    t = mdsvalue(['dim_of(\top.active_mhd.signals:BP', bp_str{j}, ',0)']);
    % time measured at 5MHz, every 0.2 us -- very fast!

    bp_dot = mdsvalue(['\top.active_mhd.signals:BP', bp_str{j}]); % [T/s]
    % also measured at 5MHz
    
    if strcmp(bp_dot(1), 'J') || strcmp(t(1), 'J') || length(bp_dot) == 1 || length(t) == 1
        % if bp_dot,t are not callable or only one point, we can't do the fit
        amp_array = zeros(1, length(timebase)) + NaN; % return NaN
        return 
    end 
    
    bp{j} = interp1(t,bp_dot,time); % only take the times we want in timebase to reduce computation time 

    for k = 1:length(time)

        b_matrix(k, j) = bp{j}(k); % rows = timebase times, columns = toroidal locations

    end 

end

mdsclose;

amplitudes = zeros(1, length(time)); % initialize

for l = 1:length(time) % for all of the rows in b_matrix = times in timebase

    M = max(b_matrix(l,:)); % find the maximum (of the toroidal amplitudes)

    x0 = [M, 0]; % use this maximum as the initial guess for the amplitude
    % of the cosine fit and 0 as the guess for the phase

    %[x,resnorm,residual,exitflag,output,lagrange,jacobian] = ... 
    x = lsqcurvefit(@n1_mode_fit_function, x0, rad, b_matrix(l,:), lb, ub, option);
    % use n1_mode_fit_function.m to fit a cosine curve to the data
    
    amplitudes(l) = x(1); % array of amplitudes for each time

end

% mdsopen('mhd', shot);
% 
% Btor = mdsvalue('-\magnetics::Btor'); % toroidal B-field
% tB = mdsvalue('dim_of(-\magnetics::Btor)'); % time
% 
% B_T = interp1(tB, Btor, time); % toroidal B-field at the times in timebase
B_T = 1;
amp_array = amplitudes./B_T; % normalize B_p,n=1 amplitudes to B_T [Hz]

if nargin == 3 % if num_times was provided, we need to find the values at the times in timebase
    
    amp_array = interp1(time, amp_array, timebase);
    
end

mdsclose;

% reorient all output arrays to original configuration
amp_array = orient_array(amp_array, original_orientation);

end