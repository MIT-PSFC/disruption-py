function [n_equal_1_amplitude_bpdot, n_equal_1_phase_bpdot] = ...
  get_n_equal_1_amplitude_bpdot(shot, timebase);

% Original author: Alex Tinguely 2015-09-14
%
% Updated 2016-03-18 by Alex Tinguely to include changes to the names of
% the the Mirnov coils, specifically KA_TOP to GH_TOP. But this would only
% apply to shots for 2016.
%
% Rewritten by R. Granetz 2016/04/28 to reformulate as a linear least
% squares fit to an expansion in terms of n = 0, 1, and 2 toroidal harmonics.
%
% Edited get_n_equal_1_amplitude.m by Alex Tinguely 2016-04-28. The script
% calculates the n=0,1,2 modes amplitudes and phases based on the time derivative
% of the poloidal magnetic field (i.e. fluctuations in time, Bp_dot [T/s]),
% instead of integrating the signal. We do this because there is a low
% frequency Bp_dot signal that causes the integrated signal to drift, and
% the drift differs for each toroidal location, leading to "false" modes
% identified.
%
% Inputs:
%   shot = shot number
%   timebase = array of requested times
%
% Outputs:
%   n_equal_1_amplitude_bpdot = array of n=1 mode amplitudes of Bp_dot vs time,
%                         (can be normalized to Btor (dBp_dt,n=1)/Btor)
%   n_equal_1_phase_bpdot     = array of n=1 phase angles vs time

% Initialize the outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.
tic;
n_equal_1_amplitude_bpdot = NaN(size(timebase));
n_equal_1_phase_bpdot     = NaN(size(timebase));

% All the signals used inside this routine are column vectors, since
% mdsvalue returns column vectors.  But the input parameter array,
% 'timebase', may be either a column or row vector, and we want to return
% the result with the same shape as the 'timebase' parameter.  Therefore,
% in order to continue with our column-based calculations, we need to
% create a copy of 'timebase' that is guaranteed to be a column vector.  So
% just before returning from this routine, the output array will be
% converted to the same shape (i.e. column or row vector) as the 'timebase'
% input parameter.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else
  timebase_column = transpose(timebase);
end

[shotopened, status] = mdsopen('magnetics', shot);
if (mod(status,2) == 0);
  return;
end;

% We use the signals from 5 toroidally-distributed Mirnov sensors in this
% routine.  The MDS nodenames of these 5 signals include the geometrical
% location of the sensors.  Between the 2015 and 2016 campaigns, one of the
% Mirnov sensors was moved to a different toroidal location, and therefore
% a new nodename was created for its signal.  Therefore the list of
% nodenames and toroidal angles for shots in 2016 is different from prior
% years.  

% ALSO, after the sensor was moved, its signal polarity was
% mistakenly flipped, and must be for data in 2016.
% NOTE by Alex Tinguely 04-28-16. This was changed by Ted G. so we comment
% out that section below.

if shot < 1160101000;  % for shots before 2016
  bp_str = {'_AB_TOP', '_BC_TOP', '_EF_TOP', '_KA_TOP', '06_GHK'}; % nodenames
  phi = [349.8; 300.1; 190.4; 15.2; 128.1];  % toroidal locations (degrees)
else;  % for shots during the 2016 campaign
  bp_str = {'_AB_TOP', '_BC_TOP', '_EF_TOP', '_GH_TOP', '06_GHK'}; 
  phi = [349.8; 300.1; 190.4; 120; 128.1];
end

phi = phi*pi/180; % change to radians

nsensors = length(phi);
ntimes = length(timebase);

bp_dot = NaN(nsensors, ntimes);

for isensor = 1:nsensors;
  nodename = ['\top.active_mhd.signals:BP' bp_str{isensor}];
  [time, status1] = mdsvalue(['dim_of(' nodename ')']);  % 5 MHz sampling!
  [Bp_dot, status2] = mdsvalue(nodename);  % tesla/s

% If the signals are not available, or only have one point, or are not the
% same length, we can't do the fit

  if mod(status1,2) == 0 || mod(status2,2) == 0 || ...
     length(Bp_dot) == 1 || length(time) == 1 || ...
     length(time) ~= length(Bp_dot);
    return; 
  end;
    
% if (shot > 1160101000 && isensor == 4);  % correct polarity for 2016 data
%   bp_dot = -bp_dot;                      % for BP_GH_TOP (the relocated coil)
% end;

%  Bp = cumtrapz(time, bp_dot);  % Integrate bp_dot [T/s] to get Bp [T]

  bp_dot(isensor,:) = interp1(time, Bp_dot, timebase_column); % Interpolate Bp_dot onto
                                                      % requested timebase 
end;

% We can normalize the Bp fluctuations by the toroidal magnetic field,
% leading to amplitudes of unit Hz. For now I will comment out this section
% to NOT normalize the data, and leave in units T/s.

% Before closing the magnetics tree, get the toroidal magnetic field.  This
% will be used to normalize the n=1 Bp component.

% btor = mdsvalue('\btor'); % tesla
% btor_time = mdsvalue('dim_of(\btor)');

mdsclose;

% btor = interp1(btor_time, btor, timebase_column); % Interpolate btor(t) onto
                                                  % the requested timebase 

% Create the 'design' matrix ('A') for the linear system of equations:
%   Bp(phi) = C1 + C2*sin(phi) + C3*cos(phi) + C4*sin(2*phi) + C5*cos(2*phi)

ncoeffs = 5;
A = NaN(nsensors, ncoeffs);

A(:,1) = [1; 1; 1; 1; 1];
A(:,2) = sin(phi);
A(:,3) = cos(phi);
A(:,4) = sin(2*phi);
A(:,5) = cos(2*phi);

% Solve the linear least squares problem to get the fitting coefficients
% for all time values
%   The n=1 amplitude at each time is sqrt(C2^2 + C3^2)
%   The n=1 phase     at each time is -atan(C3/C2)

coeffs = A \ bp_dot;
n_equal_1_amplitude_bpdot = transpose(sqrt(coeffs(2,:).^2 + coeffs(3,:).^2));
n_equal_1_phase_bpdot     = transpose(-atan(coeffs(2,:) ./ coeffs(3,:)));

% n_equal_1_amplitude = n_equal_1_amplitude ./ btor;

% If 'timebase' is a row vector, convert outputs to row vectors

if (size(timebase,2) > 1)
  n_equal_1_amplitude_bpdot   = transpose(n_equal_1_amplitude_bpdot);
  n_equal_1_phase_bpdot       = transpose(n_equal_1_phase_bpdot);
end;
toc
end