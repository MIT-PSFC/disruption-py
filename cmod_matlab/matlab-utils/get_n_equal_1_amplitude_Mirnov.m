function [n_equal_1_amplitude, n_equal_1_phase] = ...
  get_n_equal_1_amplitude_Mirnov(shot, timebase);

% With this routine I examined the high-frequency Mirnov sensors used by
% Ted Golfinopoulos to see if they could yield credible information about
% the amplitude and phase of the n = 1 component of the non-axisymmetric Bp
% field.  The signals from these sensors need to be digitially integrated,
% unlike the equilibrium BP sensors, which are analog integrated.  After
% careful study, I have concluded that the Mirnov sensors can NOT provide
% useful n = 1 information.
%
% Written by R. Granetz   2018/01
%
% Inputs:
%   shot = shot number
%   timebase = array of requested times
%
% Outputs:
%   n_equal_1_amplitude = array of n=1 mode amplitude of Bp(phi) vs time,
%                         normalized to Btor
%   n_equal_1_phase     = array of n=1 phase angles vs time

% Initialize the outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

n_equal_1_amplitude = NaN(size(timebase));
n_equal_1_phase     = NaN(size(timebase));

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

% We use the signals from 13 toroidally-distributed Mirnov sensors, all on
% the outboard side of the vessel/limiters at z = 0.108 +/- .020 m routine,
% as specified to me by Ted Golfinopoulos.  Note that the set of valid
% Mirnov signals changed between the 2015 and 2016 campaigns.

if shot < 1160101000;  % for shots before 2016
  nodename = {'BP_AB_TOP', ...
              'BP1T_ABK' , ...
              'BP2T_ABK' , ...
              'BP3T_ABK' , ...
              'BP05_ABK' , ...
              'BP_BC_TOP', ...
              'BP_EF_TOP', ...
              'BP1T_GHK' , ...
              'BP2T_GHK' , ...
              'BP3T_GHK' , ...
              'BP05_GHK' , ...
              'BP_KA_TOP'};
  phi =      [- 10.16, ...
              - 23.10, ...
              - 25.50, ...
              - 27.90, ...
              - 30.60, ...
              - 59.87, ...
              -169.55, ...
              -224.40, ...
              -226.80, ...
              -229.20, ...
              -231.90, ...
              -344.80];
  theta =    [  18.00, ...
                23.70, ...
                23.70, ...
                23.70, ...
                22.10, ...
                18.00, ...
                19.84, ...
                23.10, ...
                23.10, ...
                23.10, ...
                22.80, ...
                18.00];
else;  % for shots during the 2016 campaign
  nodename = {'BP_AB_TOP', ...
              'BP17_ABK' , ...
              'BP1T_ABK' , ...
              'BP2T_ABK' , ...
              'BP3T_ABK' , ...
              'BP05_ABK' , ...
              'BP_BC_TOP', ...
              'BP_EF_TOP', ...
              'BP20_GHK' , ...
              'BP1T_GHK' , ...
              'BP2T_GHK' , ...
              'BP3T_GHK' , ...
              'BP05_GHK' }; 
  phi =      [- 10.16, ...
              - 20.40, ...
              - 23.10, ...
              - 25.50, ...
              - 27.90, ...
              - 30.60, ...
              - 59.87, ...
              -169.55, ...
              -221.70, ...
              -224.40, ...
              -226.80, ...
              -229.20, ...
              -231.90];
  theta =    [  18.00, ...
                22.10, ...
                23.70, ...
                23.70, ...
                23.70, ...
                22.10, ...
                18.00, ...
                19.84, ...
                22.80, ...
                23.10, ...
                23.10, ...
                23.10, ...
                22.80];
end;

phi = phi*pi/180; % change to radians

nsensors = length(phi);
ntimes = length(timebase);

bp = NaN(nsensors, ntimes);

mdstcl('set default \top.active_mhd.signals');
for isensor = 1:nsensors;
  [time, status1] = mdsvalue(['dim_of(' nodename{isensor} ')']); % 5 MHz!
  [bp_dot, status2] = mdsvalue(nodename(isensor));  % tesla/s

% If the signals are not available, or only have one point, or are not the
% same length, we can't do the fit.
% 2017/12/23 -- Found a weird bug on shot 1140702005; The timebase
% of BP_BC_TOP was read in without a status error, but it was all NaN.  So
% I added some additional tests for that.

  if mod(status1,2) == 0 || mod(status2,2) == 0 || ...
     length(bp_dot) == 1 || length(time) == 1 || ...
     length(find(~isnan(time)))   <= 1 || ...
     length(find(~isnan(bp_dot))) <= 1 || ...
     length(time) ~= length(bp_dot);
    mdsclose;
    return; 
  end;
    
% I want to use the z-component of the poloidal field.  However, the Mirnov
% sensor coils are not vertical, but rather tilted poloidally by varying
% amounts, since they are mounted parallel to curved limiter shapes.  So
% they measure less than the actual Bz, and I need to compensate for that:

% Also, the polarity of the signals must be inverted to be consistent
% with +z being up.

  bp_dot = -bp_dot/cos(theta(isensor)*pi/180);

% Before integrating the signal, subtract any DC offset due to the
% digitizer.  This is necessary to prevent the integrated signal from
% effectively having a baseline drift.  Ted instructed me to use all times
% before t = 0 s to determine the baseline.

  indices_baseline = find(time < 0);
  baseline = mean(bp_dot(indices_baseline));
  bp_dot = bp_dot - baseline;

% Now integrate bp_dot [T/s] to get Bp [T]

  Bp = cumtrapz(time, bp_dot);

% Interpolate Bp onto requested timebase 

  bp(isensor,:) = interp1(time, Bp, timebase_column);
end;

% Before closing the magnetics tree, get the toroidal magnetic field.  This
% will be used to normalize the n=1 Bp component.

btor_time = mdsvalue('dim_of(\btor)');
btor = mdsvalue('\btor'); % tesla
polarity = sign(sum(btor));
if polarity == -1;
  btor = -btor;
end;

mdsclose;

% Interpolate btor(t) onto the requested timebase 

btor = interp1(btor_time, btor, timebase_column);

%{

% Create the 'design' matrix ('A') for the linear system of equations:
%   Bp(phi,t) = C1(t) + C2(t)*sin(phi) + C3(t)*cos(phi) + ...
%               C4(t)*sin(2*phi) + C5(t)*cos(2*phi)

ncoeffs = 5;
A = NaN(nsensors, ncoeffs);

A(:,1) = ones(nsensors, 1);
A(:,2) = sin(phi);
A(:,3) = cos(phi);
A(:,4) = sin(2*phi);
A(:,5) = cos(2*phi);

% Solve the linear least squares problem to get the fitting coefficients
% for all time values

coeffs = A \ bp;

%}

% Create the 'design' matrix ('A') for the linear system of equations:
%   Bp(phi,t) = C1(t) + C2(t)*sin(phi) + C3(t)*cos(phi)

ncoeffs = 3;
A = NaN(nsensors, ncoeffs);

A(:,1) = ones(nsensors, 1);
A(:,2) = sin(phi);
A(:,3) = cos(phi);

% Solve the linear least squares problem to get the fitting coefficients
% for all time values

coeffs = A \ bp;

% The n=1 amplitude at each time is sqrt(C2^2 + C3^2)
% The n=1 phase     at each time is -atan(C2/C3)

n_equal_1_amplitude = transpose(sqrt(coeffs(2,:).^2 + coeffs(3,:).^2));
n_equal_1_phase     = transpose(-atan(coeffs(2,:) ./ coeffs(3,:)));

% Normalize the n=1 amplitude to Btor

n_equal_1_amplitude = n_equal_1_amplitude ./ btor;

% If 'timebase' is a row vector, convert outputs to row vectors

if (size(timebase,2) > 1)
  n_equal_1_amplitude   = transpose(n_equal_1_amplitude);
  n_equal_1_phase       = transpose(n_equal_1_phase);
end;

% For debugging purposes only, use the fitted coefficients to reconstruct
% the signals from the n=0, and 1 components.  This section will only be
% done if test=1.  I can do this by running this routine with the
% debugger, stopping at the following statement, and manually setting
% test=1.  Otherwise this section will be skipped.

if exist('test','var') && ~isempty(test);
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
