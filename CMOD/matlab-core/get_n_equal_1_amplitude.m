function [n_equal_1_amplitude, n_equal_1_normalized, n_equal_1_phase] = ...
  get_n_equal_1_amplitude(shot, timebase);

% Calculate n=1 amplitude and phase for the disruption warning database
% using the four BP13 Bp sensors near the midplane on the outboard vessel
% wall.  The calculation is done by using a least squares fit to an
% expansion in terms of n = 0 & 1 toroidal harmonics.  The BP13 sensors are
% part of the set used for plasma control and equilibrium reconstruction,
% and their signals have been analog integrated (units: tesla), so they
% don't have to be numerically integrated.  These four sensors were working
% well in 2014, 2015, and 2016.  I looked at our locked mode MGI run on
% 1150605, and the different applied A-coil phasings do indeed show up on
% the n=1 signal.
%
% I have also examined the set of outboard Mirnov coils that Ted
% Golfinopoulos is using.  These have to be numerically integrated.  I
% found that there is too much variation in their signal amplitudes to
% allow for a credible extraction of n = 1 information.
% 
% Written by R. Granetz 2016/05/04 
% Revision history:
%    2018/01    -- R. Granetz; some modifications made
%
% Inputs:
%   shot = shot number
%   timebase = array of requested times
%
% Outputs:
%   n_equal_1_amplitude  = array of n=1 mode amplitude of Bp vs time [tesla]
%   n_equal_1_normalized = array of n=1 mode amplitude of Bp vs time,
%                            normalized to Btor
%   n_equal_1_phase      = array of n=1 phase angles vs time [radians]

% Initialize the outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

n_equal_1_amplitude  = NaN(size(timebase));
n_equal_1_normalized = NaN(size(timebase));
n_equal_1_phase      = NaN(size(timebase));

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

% We use the signals from 4 toroidally-distributed BP13 sensors, all on
% the outboard side of the vessel wall at R = 1.0573 m, z = 0.0845 m.
% But their toroidal field pickup must be subtracted and their baselines
% need to be subtracted.

bp13_names = {'BP13_BC', 'BP13_DE', 'BP13_GH', 'BP13_JK'};
nsensors = length(bp13_names);
ntimes = length(timebase);
bp13_signals = NaN(nsensors, ntimes);

mdstcl('set default \mag_bp_coils');
bp_nodenames = mdsvalue('nodename');
phi = mdsvalue('phi');
btor_pickup_coeffs = mdsvalue('btor_pickup');

[~,bp13_indices,~] = intersect(bp_nodenames, bp13_names);

bp13_phi = phi(bp13_indices) + 360;
bp13_btor_pickup_coeffs = btor_pickup_coeffs(bp13_indices);

btor = mdsvalue('\btor');
magtime = mdsvalue('dim_of(\btor)');
baseline_indices = find(magtime <= -1.8);
btor = btor - mean(btor(baseline_indices));

mdstcl('set default \mag_bp_coils.signals');

for isensor = 1:nsensors;
  [signal, status] = mdsvalue(bp13_names(isensor));  % tesla

% If the signals are not available, or only have one point, or are not the
% same length as their timebase, we can't do the fit

  if mod(status,2) == 0 || length(signal) == 1 || ...
      length(magtime) ~= length(signal);
    return; 
  end;
    
% Subtract baseline offset

  baseline = mean(signal(baseline_indices));
  signal = signal - baseline;

% Subtract Btor pickup

  signal = signal - bp13_btor_pickup_coeffs(isensor) * btor;

% Interpolate Bp onto requested timebase 

  bp13_signals(isensor,:) = interp1(magtime, signal, timebase_column);
end;
mdsclose;

% Interpolate btor_magnitude(t) onto the requested timebase 

polarity = sign(mean(btor));
btor_magnitude = btor * polarity;
btor_magnitude = interp1(magtime, btor_magnitude, timebase_column);

% Create the 'design' matrix ('A') for the linear system of equations:
%   Bp(phi) = C1 + C2*sin(phi) + C3*cos(phi)

ncoeffs = 3;
A = NaN(nsensors, ncoeffs);

A(:,1) = transpose([1,1,1,1]);
A(:,2) = transpose(sin(bp13_phi*pi/180));
A(:,3) = transpose(cos(bp13_phi*pi/180));

% Solve the linear least squares problem to get the fitting coefficients
% for all time values

coeffs = A \ bp13_signals;

% The n=1 amplitude at each time is sqrt(C2^2 + C3^2)
% The n=1 phase     at each time is arctan(-C2/C3), using complex number
%   phasor formalism, exp(i(phi - delta))

n_equal_1_amplitude = transpose(sqrt(coeffs(2,:).^2 + coeffs(3,:).^2));
n_equal_1_phase     = transpose(atan2(-coeffs(2,:), coeffs(3,:)));

% Normalize the n=1 amplitude to Btor

n_equal_1_normalized = n_equal_1_amplitude ./ btor_magnitude;

% If 'timebase' is a row vector, convert outputs to row vectors

if (size(timebase,2) > 1)
  n_equal_1_amplitude   = transpose(n_equal_1_amplitude);
  n_equal_1_normalized  = transpose(n_equal_1_normalized);
  n_equal_1_phase       = transpose(n_equal_1_phase);
end;

% For debugging purposes only, use the fitted coefficients to reconstruct
% the signals from the n=0, and 1 components.  This section will only be
% done if test=1.  I can do this by running this routine with the
% debugger, stopping at the following statement, and manually setting
% test=1.  Otherwise this section will be skipped.

if exist('test','var') && test == 1;
  phi_array = [0 : 1 : 360] * pi/180;

  A_reconstruct = NaN(length(phi_array), ncoeffs);

  A_reconstruct(:,1) = ones(length(phi_array), 1);
  A_reconstruct(:,2) = sin(phi_array);
  A_reconstruct(:,3) = cos(phi_array);
% A_reconstruct(:,4) = sin(2*phi_array);
% A_reconstruct(:,5) = cos(2*phi_array);

  bp_reconstructed = A_reconstruct * coeffs;  
end;

end
