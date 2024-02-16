function v_loop = get_v_loop(shot, timebase);

% This routine gets the loop voltage from the EAST tree.  The signal in the
% tree is derived by taking the time derivative of a flux loop near the
% inboard midplane.  Two possible signals are available, one digitised at a
% high rate (50 kHz), and the other sub-sampled down to 1 kHz.  This
% routine reads in the 1 kHz signal.  It linearly interpolates the loop
% voltage signal onto the specified timebase.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Output:
%   v_loop = loop voltage [V]
%
% Author: Robert Granetz   Apr 2016

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

% Initialize all output arrays to NaN (Not-a-Number) column vectors

v_loop = NaN(length(timebase), 1);

% Get "\vp1_s" signal from the EAST tree.  (This signal is a sub-sampled
% version of "vp1".)

[shotopened, status] = mdsopen('east', double(shot));
if (mod(status,2) == 1);
  [v_loop, status] = mdsvalue('\vp1_s');  % Read in loop voltage data [volts]
  if (mod(status,2) == 1);                % If successful, continue
    v_loop_time = mdsvalue('dim_of(\vp1_s)');
% Interpolate the signal onto the requested timebase.
    v_loop = interp1(v_loop_time, v_loop, timebase_column, 'linear');
  else;
    v_loop = NaN(length(timebase), 1);
  end;
  mdsclose;
end;

% The output signal is currently a column vector.  However, we desire to
% have the output array match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert the output to a row
% vector.

if (size(timebase,2) > 1);
  v_loop = transpose(v_loop);
end;

end
