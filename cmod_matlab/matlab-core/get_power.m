function [p_RAD, dprad_dt, p_LH, p_OHM, p_ICRF, p_input, rad_fraction, ...
  v_loop] = get_power(shot, timebase);

% This function gets the input heating powers -- ohmic (p_OHM), ion
% cyclotron (p_ICRF), and lower hybrid (p_LH) -- as well as the radiated
% output power (p_RAD).  If any of the auxiliary heating powers are not
% available (there was no ICRF or LH), then this function returns an array
% of zeros for them.  The ohmic heating power is obtained by calling a
% separate function, get_P_ohm.m.  If either the ohmic heating power or the
% radiated power is unavailable, then arrays of NaN (Not-a-Number) are
% returned for them, and for the radiated power fraction.
%
% Revision/history:
%  2016/02    - RSG; Changed p_RAD to use the \twopi_diode signal rather
%                     than \twopi_foil.  This is because the processing of
%                     the foil bolometer signal involved an unacceptable
%                     degree of non-causal filtering.  The diode signal has
%                     been cross-calibrated against the foil bolometer
%                     during flattop periods of many non-disruptive
%                     discharges.
%  2016/03/15 - RSG; Added dprad_dt

% Inputs:
%   shot = shot number
%   timebase = times to calculate the powers 
% Outputs:
%   p_RAD  = radiated power [W]
%   dprad_dt = time derivative of p_RAD [W/s]
%   p_LH   = lower hybrid power [W]
%   p_OHM  = ohmic power [W]
%   p_ICRF = ion cyclotron power [W]
%   p_input = total input power = p_OHM + p_LH + p_ICRF [W]
%   rad_fraction = p_RAD/(p_OHM + p_LH + p_ICRF), ratio of radiated output
%     power to input heating power 
%   v_loop = edge loop voltage; time derivative of a weighted average of
%     flux loops obtained from MFLUX
%
% Author: Alex Tinguely 15-09-16
% Some modifications by Robert Granetz 2015/10/11

% We want the output vectors to have the same shape as the input vector,
% "timebase", i.e. if "timebase" is a row vector, the output vectors will
% be row vectors; if "timebase" is a column vector, the output vectors will
% be column vectors.  Note that Matlab calls to mdsvalue return column
% vectors for 1-D signals, so this routine works internally with column
% vectors.  So the first step is to create a column vector copy of the
% input parameter, "timebase".

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% ------------------------------------
% Get lower hybrid power
%
% NOTE: the timebase for the LH power signal does not extend over the full
% time span of the discharge.  Therefore, when interpolating the LH power
% signal onto the "timebase" array, the LH signal has to be extrapolated
% with zero values.  This is an option in the 'interp1' routine.  If the
% extrapolation is not done, then the 'interp1' routine will assign NaN
% (Not-a-Number) values for times outside the LH timebase, and the NaN's
% will propagate into p_input and rad_fraction, which is not desirable.

[shotopened, status] = mdsopen('LH', shot);

if (mod(status,2)==0);
  p_LH = zeros(size(timebase));
else;
  p_lh = 1.e3 * mdsvalue('\LH::TOP.RESULTS:NETPOW'); % [W]
  [t_lh, status] = mdsvalue('dim_of(\LH::TOP.RESULTS:NETPOW)'); % [s]
  if (mod(status,2)==0 || length(t_lh) == 1);
    p_LH = zeros(size(timebase));
  else;
    p_LH = interp1(t_lh, p_lh, timebase_column, 'linear', 0.);
    if (size(timebase,2) > 1); % If "timebase" is a row vector, then make
      p_LH = transpose(p_LH);  % p_LH into a row vector.
    end;
  end;
  mdsclose;
end;

% ------------------------------------
% Get ICRF power

[shotopened, status] = mdsopen('RF', shot);

if (mod(status,2)==0);
  p_ICRF = zeros(size(timebase));
else;
  p_icrf = 1.e6 * mdsvalue('\rf::rf_power_net'); % [W]
  [t_icrf, status] = mdsvalue('dim_of(\rf::rf_power_net)'); % [s]
  if (mod(status,2)==0 || length(t_icrf) == 1);
    p_ICRF = zeros(size(timebase));
  else;
    p_ICRF = interp1(t_icrf, p_icrf, timebase_column, 'linear', 0.);
    if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
      p_ICRF = transpose(p_ICRF); % p_ICRF into a row vector.
    end;
  end;
  mdsclose;
end;

% ------------------------------------
% Get ohmic power and loop voltage

[p_OHM, v_loop] = get_P_ohm(shot, timebase); %[W, V]

% ------------------------------------
% Radiated power

[shotopened, status] = mdsopen('spectroscopy', shot);
if (mod(status,2)==0);
  p_rad = NaN(size(timebase));
  dprad_dt = NaN(size(timebase));
else;
% p_rad = mdsvalue('\twopi_foil'); % [W]
% [t_rad, status] = mdsvalue('dim_of(\twopi_foil)'); % [s]
  p_rad = mdsvalue('\twopi_diode') * 1.e3; % [W]
  p_rad = p_rad * 4.5; % Factor of 4.5 comes from cross-calibration with
                       % 2pi_foil during flattop times of non-disruptive
                       % shots, excluding times for
                       % which p_rad (uncalibrated) <= 1.e5 W
  [t_rad, status] = mdsvalue('dim_of(\twopi_diode)'); % [s]
  if (mod(status,2)==0 || length(t_rad) == 1);
    p_RAD = NaN(size(timebase));
    dprad_dt = NaN(size(timebase));
  else;   
    p_RAD = interp1(t_rad, p_rad, timebase_column, 'linear');
    dprad_dt = gradient(p_rad, t_rad); % No smoothing for now
    dprad_dt = interp1(t_rad, dprad_dt, timebase_column, 'linear');
    if (size(timebase,2) > 1);  % If "timebase" is a row vector, then make
      p_RAD = transpose(p_RAD); % p_RAD into a row vector.
      dprad_dt = transpose(dprad_dt);
    end;
  end;
  mdsclose;
end

% ------------------------------------
p_input = p_OHM + p_LH + p_ICRF; % Total input power [W]
rad_fraction = p_RAD./p_input;   % Ratio of output power to input power

rad_fraction(find(isinf(rad_fraction))) = NaN; % Change any infinite values
                                               % to Not-a-Number
end
