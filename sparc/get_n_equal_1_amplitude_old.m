function [n_equal_1_amp] = get_n_equal_1_amplitude_old(shot, timebase, num_times)

% This function finds the best fit of the n = 1  MHD mode around the torus 
% for five different toroidal locations (Mirnov coils) using fluctuations 
% in the TIME DERIVATIVE of the poloidal magnetic field (T/s). This is  
% based off of work done by Geoff Olynyk in his thesis. Then this is  
% normalized by the toroidal magnetic field, so the resulting n=1 amplitude 
% ratio is unitless.
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
% Updated 2016-03-18 by Alex Tinguely to include changes to the names of
% the the Mirnov coils, specifically KA_TOP to GH_TOP. But this would only
% apply to shots for 2016.
%
% Inputs:
%   shot = shot number
%   timebase = array of requested times
%   num_times = OPTIONAL ARGUMENT, number of time slices to evaluate the 
%      n=1 mode fit. A greater number of divisions (large num_times) means
%               a better fit. If no argument is provided, it will default
%               to the number of elements in timebase.
%
% Outputs:
%   n_equal_1_amp = array of n=1 mode amplitudes in B_poloidal =
%   (dB_p,n=1)/B_T
%tic;
format long % increase the precision of the numbers in the mdsvalues

% Initialize the output to NaN (Not-a-Number).  This will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

n_equal_1_amp = NaN(size(timebase));

% If the third argument, num_times, is provided, we will 

if nargin == 3
    time = linspace(timebase(1), timebase(end), num_times);    
else 
    time = timebase; 
end

% All the signals used inside this routine are column vectors, since
% mdsvalue returns column vectors.  But the input parameter, 'timebase',
% may be either a column or row vector.  Therefore, in order to continue
% with our column-based calculations, we need to create a copy of
% 'timebase' that is guaranteed to be a column vector.  Just before
% returning from this routine, all output vectors will be converted to the
% same shape (i.e. column or row vector) as the 'timebase' input
% parameter.

if (size(time,1) > 1);
  timebase_column = time;
else
  timebase_column = transpose(time);
end

[shotopened, status] = mdsopen('magnetics', shot);
if (mod(status,2) == 0);
    return;
end;

% Before the 2016 campaign, we used five Mirnov coils with names listed in
% the first if statement, but _KA_TOP has since become _GH_TOP, so we need
% to grab the appropriate data depending on the year

if shot < 1160101000 % for shots before 2016
    bp_str = {'_AB_TOP', '_BC_TOP', '_EF_TOP', '_KA_TOP', '06_GHK'}; 
    % 5 Mirnov coil names
    phi = [349.8, 300.1, 190.4, 15.2, 128.1]; 
    % 5 toroidal locations (degrees)
else % for shots during the 2016 campaign
    bp_str = {'_AB_TOP', '_BC_TOP', '_EF_TOP', '_GH_TOP', '06_GHK'}; 
    phi = [349.8, 300.1, 190.4, 120, 128.1];
end

phi = phi*pi/180; % change to radians

% Now we set the bounds for the fit function below.
lb = [0, -pi]; % lower bound, don't want negative amplitudes
ub = [Inf, pi]; % upper bound 

% In the options, increase the number of function evaluations and suppress
% the display of text after every iteration
option = optimset('TolFun', 0, 'MaxFunEvals', 1e8, 'Display', 'off');

bp = cell(1,5); % initialize cell array for all 5 locations
b_matrix = zeros(length(timebase_column), 5);

for i = 1:length(bp_str) % for each toroidal location

    [t, status] = mdsvalue(['dim_of(\top.active_mhd.signals:BP',...
        bp_str{i}, ',0)']);
    % time measured at 5MHz, every 0.2 us -- very fast!

    [bp_dot, status] = mdsvalue(['\top.active_mhd.signals:BP',...
        bp_str{i}]); % [T/s], also measured at 5MHz
    
    % if bp_dot,t are not callable, only one point, 
    % or not the same length, we can't do the fit
    if mod(status,2) == 0 || length(bp_dot) == 1 || length(t) == 1 || ...
            length(t) ~= length(bp_dot)
        return; 
    end
    
    %size(t)
    %size(bp_dot)
    
    % cumulative numerical integration (trapezoidal) [T]
    bp1 = cumtrapz(t,bp_dot);
    
    bp{i} = interp1(t,bp1,timebase_column); 
    % only take the times we want in timebase to reduce computation time 

    for j = 1:length(timebase_column)
        b_matrix(j, i) = bp{i}(j); 
        % rows j = timebase times, columns i = toroidal locations
    end 
end

mdsclose;

amplitudes = zeros(length(timebase_column),1); % initialize
%toc;
% for all of the rows in b_matrix = times in timebase
%tic;
for i = 1:length(timebase_column) 

    M = max(b_matrix(i,:)); % find the maximum (of the toroidal amplitudes)

    x0 = [M, 0]; % use this maximum as the initial guess for the amplitude
    % of the cosine fit and 0 as the guess for the phase

    x = lsqcurvefit(@n1_mode_fit_function, x0, phi, b_matrix(i,:), ...
        lb, ub, option);
    % use n1_mode_fit_function.m to fit a cosine curve to the data
    
    amplitudes(i) = x(1); % array of amplitudes for each time
end

% Now we get the toroidal magnetic field.
[shotopened, status] = mdsopen('mhd', shot);
if (mod(status,2) == 0);
    return;
end;

Btor = mdsvalue('-\magnetics::Btor'); % toroidal B-field
tB = mdsvalue('dim_of(-\magnetics::Btor)'); % time

B_T = interp1(tB, Btor, timebase_column); 
% toroidal B-field at the times in timebase

n_equal_1_amp = amplitudes./B_T; % normalize B_p,n=1 amplitudes to B_T

% if num_times was provided, we need to find the values at the times in 
% timebase
if nargin == 3 
    n_equal_1_amp = interp1(timebase_column, n_equal_1_amp, timebase);
end

mdsclose;

% If timebase is a row vector, output n_equal_1_amp as a row vector.
if (size(timebase,2) > 1)
    n_equal_1_amp = transpose(n_equal_1_amp);
end
%toc;
end