function [amplitude, phase, time] = ...
  get_saddle_fourier_components(shot, do_not_subtract_RMP);

if ~exist('do_not_subtract_RMP', 'var');
  do_not_subtract_RMP = 0;
elseif isempty(do_not_subtract_RMP);
  do_not_subtract_RMP = 0;
elseif do_not_subtract_RMP ~= 1;
  do_not_subtract_RMP = 0;
end;

% This routine reads in the data from the toroidal array of 8 saddle coils
% on EAST, and then does a fast Fourier transform (FFT) to determine the
% amplitude and phase of the n = 0, 1, 2, and 3 components as a function
% of time.
%
%  2016/12/05 -- Gu Xiang: write Matlab code to read in saddle coil signals
%                  and calculate Fourier component amplitudes.
%  2016/12/07 -- R. Granetz: modify Gu Xiang's code to convert it to a
%                  function, add the phase information, improve
%                  computational efficiency, and add comments
%

% -------------------------------------------------------------------
%   Check the lock mode
%   05/12/2016 by Xiang Gu
% -------------------------------------------------------------------

% Read in the saddle sensor data and the rmp currents from the MDSplus
% tree.  All the outputs have time as their 1st dimension,
% i.e. rmp(time,rmpcoil#) and saddle(time,saddlecoil#), and the time arrays
% are column vectors.

% The ordering of the rmp coils is: irmpu1 to irmpu8, then irmpl1 to irmpl8.
% The ordering of the saddle coils is: sad_pa, sad_bc, .. sad_no.

[rmptime, rmp, saddletime, saddle] = get_rmp_and_saddle_signals(shot);

if length(saddletime) <= 2;
  amplitude = [];
  phase = [];
  time = [];
  return;
end;

if (do_not_subtract_RMP ~= 1);

% The saddle signals can include direct pickup from the RMP coils, when the
% RMP coils are active.  This direct pickup can be subtracted from the
% saddle signals, leaving just the plasma response and baseline drifts.
% The signal on each saddle sensor that is due to each RMP coil was derived
% from two RMP calibration shots, 61846 and 62053, using the Matlab routine
% "calculate_rmp_saddle_coeff_matrix.m"

  load('rmp_saddle_coeff_matrix.mat');

% Calculate the direct RMP contribution to the saddle sensors, and subtract
% it from the saddle signals.

  rmp_pickup = transpose(coeff_matrix * transpose(rmp));
  rmp_pickup = interp1(rmptime, rmp_pickup, saddletime);
  saddle = saddle - rmp_pickup;
end;

% Calculate fast Fourier transforms and get mode amplitudes and phases

output = fft(saddle, [], 2);  % Take FFT along 2nd dimension (phi)

amplitude = abs(output)/size(output,2);
amplitude(:, 2:end) = amplitude(:, 2:end) * 2;

phase = atan2(imag(output), real(output));

% Only want n = 0, 1, 2, 3 (first four Fourier modes)

amplitude = amplitude(:,1:4);
phase = phase(:, 1:4);
time = saddletime;

return;
end
