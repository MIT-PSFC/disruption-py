function [n_equal_1_mode, n_equal_2_mode, n_equal_3_mode] = ...
  get_locked_mode_data_EAST(shot, timebase);

% This routine gets the n=1, 2, and 3 locked mode amplitudes from the EAST
% saddle coil measurements.  If the error field correction coils are on,
% they can contribute directly to the saddle coil signals, so maybe we
% should subtract out their contribution.  But I have not done that yet.
% 
% Author: Robert Granetz   Dec 2016

% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, the "mdsvalue" routine in Matlab returns
% column vectors for 1-D signals, so it is simpler to work with column
% vectors within this routine, and then, if "timebase" is a row vector,
% convert the outputs to row vectors just before exiting this routine.  So
% the first step is to create a copy of "timebase" that is guaranteed to be
% a column vector.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% EAST has a toroidal array of 8 saddle coil sensors.  Gu Xiang and I have
% written a routine which gets the saddle coil signals from MDSplus, does a
% Fourier decompostion, and returns the n=0, 1, 2, & 3 amplitudes and
% phases.
%
% Note: the variables returned from our routine are:
%   amplitude(time, mode) % mode = 1 to 4, which corresponds to n = 0 to 3
%   phase(time, mode) %  -pi to +pi
%   time % This is the timebase of the saddle coil signals and the Fourier data

[amplitude, phase, time] = get_saddle_fourier_components(shot);

if (length(time) <= 1 || isnan(time(1)));
  n_equal_1_mode = NaN(length(timebase), 1);
  n_equal_2_mode = NaN(length(timebase), 1);
  n_equal_3_mode = NaN(length(timebase), 1);
  return;
end;

% Interpolate mode amplitude onto the requested timebase

n_equal_1_mode = interp1(time, amplitude(:,2), timebase_column, 'linear');
n_equal_2_mode = interp1(time, amplitude(:,3), timebase_column, 'linear');
n_equal_3_mode = interp1(time, amplitude(:,4), timebase_column, 'linear');

% The output signals are currently column vectors.  However, we desire to
% have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert the outputs to row
% vectors.

if (size(timebase,2) > 1);
  n_equal_1_mode = transpose(n_equal_1_mode);
  n_equal_2_mode = transpose(n_equal_2_mode);
  n_equal_3_mode = transpose(n_equal_3_mode);
end;

end
