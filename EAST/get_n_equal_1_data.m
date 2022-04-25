function [n_equal_1_normalized, n_equal_1_mode, n_equal_1_phase, ...
  rmp_n_equal_1, rmp_n_equal_1_phase, btor] = ...
  get_n_equal_1_data(shot, timebase);

% This function computes the amplitude and phase of the n=1 Fourier
% component of the net saddle signals (total saddle signals minus the
% calculated pickup from the RMP coils) and interpolates the signals onto
% the specified timebase.  It also computes the amplitude and phase of the
% n=1 Fourier component of the calculated pickup from the RMP coils into
% the database.
% 
% Inputs:
%   shot = shot number
%   timebase = times to sample the data at
% Outputs:
%   n_equal_1_normalized = 'n_equal_1_mode' normalized to Btor
%   n_equal_1_mode = amplitude of n=1 Fourier component of saddle signals
%                    after subtracting the calculated pickup from the RMP
%                    coil currents [tesla]
%   n_equal_1_phase = toroidal phase angle of above [radians]
%   rmp_n_equal_1 = amplitude of the n=1 Fourier component of the
%                   calculated pickup of the RMP coils on the saddle
%                   signals [tesla]
%   rmp_n_equal_1_phase = toroidal phase angle of above [radians]
%
% Author: Robert Granetz   2017/04/05
% Revision/history:
%  2018/10/11 - Chen Dalong, K. Montes, RSG: added n_equal_1_normalized and
%               associated code for Bt (Bt=\mu*(130turns*16)*it/2*\pi*R)
%  2019/11/06 - RSG; Corrected the section that gets Btor (for
%               normalization).  Starting with shot 65573 the \it signal
%               was switched back to the eng_tree.

% Set default values

n_equal_1_normalized = NaN(size(timebase));
n_equal_1_mode = NaN(size(timebase));
n_equal_1_phase = NaN(size(timebase));
rmp_n_equal_1 = NaN(size(timebase));
rmp_n_equal_1_phase = NaN(size(timebase));
btor = NaN(size(timebase));

% Read in the saddle sensor data and the rmp coil currents

[rmptime, rmp, saddletime, saddle] = get_rmp_and_saddle_signals(shot);

% If no valid saddle data, bail out

if length(saddletime) <= 2;
  return;
end;

% Calculate RMP n=1 Fourier component amplitude and phase (on the timebase
% of the saddle signals)

if length(rmptime) > 2;

% Calculate the signal on each saddle sensor due to each RMP coil current.
% The coefficient matrix that is used in this computation was derived from
% two RMP no-plasma calibration shots, 61846 and 62053, using the Matlab
% routine "calculate_rmp_saddle_coeff_matrix.m".  It is stored in the
% file 'rmp_saddle_coeff_matrix.mat'.  Since this routine is meant to be
% called many, many times to process a long list of shots, it is most
% efficient if the matrix is only read in once, and then saved between
% invocations of this routine.  This is accomplished by using the
% 'persistent' property, and testing to see if the matrix is empty or not.

  persistent coeff_matrix

  if isempty(coeff_matrix);
    load('rmp_saddle_coeff_matrix.mat');
  end;

  rmp_pickup = transpose(coeff_matrix * transpose(rmp));

% Interpolate onto the timebase of the saddle signals

  rmp_pickup = interp1(rmptime, rmp_pickup, saddletime);

% Calculate fast Fourier transforms of the RMP-induced signals, and get
% mode amplitudes and phases

  output = fft(rmp_pickup, [], 2);  % Take FFT along 2nd dimension (phi)

  amplitude = abs(output)/size(output,2);
  amplitude(:, 2:end) = amplitude(:, 2:end) * 2;

  phase = atan2(imag(output), real(output));

% Only want n = 1 Fourier component

  rmp_n_equal_1 = amplitude(:,2);
  rmp_n_equal_1_phase = phase(:,2);

% Force the RMP signals to be the same shape as the specified timebase
% (i.e. column vector or row vector)

  if (size(timebase, 2) > 1);
    saddletime = transpose(saddletime);
    rmp_n_equal_1 = transpose(rmp_n_equal_1);
    rmp_n_equal_1_phase = transpose(rmp_n_equal_1_phase);
  end;

% Interpolate onto the specified timebase

  rmp_n_equal_1 = interp1(saddletime, rmp_n_equal_1, timebase);
  rmp_n_equal_1_phase = interp1(saddletime, rmp_n_equal_1_phase, timebase);
else;
  rmp_pickup = zeros(size(saddle));
  rmp_n_equal_1 = zeros(size(timebase));
  rmp_n_equal_1_phase = NaN(size(timebase));
end;

if length(saddletime) > 2;

% The saddle signals can include direct pickup from the RMP coils, when the
% RMP coils are active.  This direct pickup must be subtracted from the
% saddle signals to leave just the plasma response and baseline drifts.

  saddle = saddle - rmp_pickup;

% Calculate fast Fourier transforms and get mode amplitudes and phases

  output = fft(saddle, [], 2);  % Take FFT along 2nd dimension (phi)

  amplitude = abs(output)/size(output,2);
  amplitude(:, 2:end) = amplitude(:, 2:end) * 2;

  phase = atan2(imag(output), real(output));

% Only want n = 1 Fourier component

  n_equal_1_mode = amplitude(:,2);
  n_equal_1_phase = phase(:, 2);

% Force the saddle signals to be the same shape as the specified timebase
% (i.e. column vector or row vector)

  if (size(timebase, 2) > 1);
    n_equal_1_mode = transpose(n_equal_1_mode);
    n_equal_1_phase = transpose(n_equal_1_phase);
  end;

% Interpolate onto the specified timebase

  n_equal_1_mode = interp1(saddletime, n_equal_1_mode, timebase);
  n_equal_1_phase = interp1(saddletime, n_equal_1_phase, timebase);

else;
  n_equal_1_mode = NaN(size(timebase));
  n_equal_1_phase = NaN(size(timebase));
end;

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% Next, get the toroidal magnetic field.  For shots < 60000, the TF current
% was in node "\it" in the "eng_tree", and the timebase of the signal was
% not well-defined.  For shots between 60000 (on 2016/01/29) and 65165
% (2016/05/01), the "\it" node is in the "pcs_east" tree, with a proper
% timebase.  For shots after 65165, the "\it" node was back in the
% "eng_tree" tree (with a proper timebase).  Also, starting with shot
% 60000, there is a fibreoptic-based measurement of the TF current.  The
% signals are "\focs_it" (digitized at 50 kHz) and "\focs_it_s"
% (sub-sampled at 1 kHz), both in the "east" tree.  The latter two signals
% differ by 1.6% from the \it signal (as of 2016/04/18).

if (shot > 65165);
  [~, status] = mdsopen('eng_tree', double(shot));
  if (mod(status, 2) == 1);
    [btor_time, status] = mdsvalue('dim_of(\it)');
    if (mod(status, 2) == 1);
      itf = mdsvalue('\it'); % in amps
      btor = (4*pi*1e-7) * itf *(16*130) /(2*pi*1.8); % about 4327 amps/tesla
      btor = interp1(btor_time, btor, timebase_column);
    else;
      btor = NaN(length(timebase), 1);
    end;
    mdsclose;
  else;
    btor = NaN(length(timebase), 1);
  end;
elseif (shot >= 60000);
  [~, status] = mdsopen('pcs_east', double(shot));
  if (mod(status,2) == 1);
    [btor_time, status] = mdsvalue('dim_of(\it)');
    if (mod(status, 2) == 1);
      itf = mdsvalue('\it'); % in amps
      btor = (4*pi*1e-7) * itf *(16*130) /(2*pi*1.8); % about 4327 amps/tesla
      btor = interp1(btor_time, btor, timebase_column);
    else;
      btor = NaN(length(timebase), 1);
    end;
    mdsclose;
  else;
    btor = NaN(length(timebase), 1);
  end;
else;   % get itf from eng_tree for shots < 60000
  [~, status] = mdsopen('eng_tree', double(shot));
  if (mod(status, 2) == 1);
    [itf, status] = mdsvalue('\it'); % in amps
    if (mod(status, 2) == 1);
      itf = mean(itf);  % Btor is constant in time (superconducting magnet)
      btor = (4*pi*1e-7) * itf *(16*130) /(2*pi*1.8); % about 4327 amps/tesla
      btor = [btor; btor];   % construct 2-point signal from scalar value
      btor_time = [0; 1000]; % construct 2-point timebase
      btor = interp1(btor_time, btor, timebase_column);
    else;
      btor = NaN(length(timebase), 1);
    end;
    mdsclose;
  else;
    btor = NaN(length(timebase), 1);
  end;
end;

btor = transpose(btor);
n_equal_1_normalized = n_equal_1_mode ./ abs(btor);

end
