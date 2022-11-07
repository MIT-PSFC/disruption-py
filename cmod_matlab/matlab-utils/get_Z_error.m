function [Z_error, Z_prog, Z_cur] = get_Z_error(shot, timebase)

% This script calculates Z_error = Z_programmed - Z_cur, i.e. how much the
% actual vertical position differs from the requested position.  It
% linearly interpolates both the programmed and measured positions onto the
% given timebase, and then takes the difference.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   Z_error = Z_programmed - Z_cur [m]
%   Z_prog = programmed (i.e. requested) vertical position of the current
%     centroid [m]
%   Z_cur = actual vertical position of the current centroid, as calculated
%     by EFIT [m]
%
% Author: Alex Tinguely 2015-09-09
% Rewritten by Robert Granetz 2016/01/29

% \hybrid::top.hardware.dpcs.signals:a_out
% \hybrid::top.dpcs.seg_??:p_02:predictor:factor  (for ZCUR)
%           or      seg_??:p_16:predictor:factor  (for IP)

% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, Matlab calls to the routine "mdsvalue"
% return column vectors for 1-D signals, so it is simpler to work with
% column vectors within this routine, and then, if necessary, convert the
% outputs to row vectors just before exiting this routine.  So the first
% step is to create a copy of "timebase" that is guaranteed to be a column
% vector.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% Initialize all output arrays to NaN (Not-a-Number) column vectors

Z_error = NaN(length(timebase),1);
Z_prog  = NaN(length(timebase),1);
Z_cur   = NaN(length(timebase),1);

% Read in the actual Z_cur, as calculated by EFIT, interpolated onto the
% requested timebase.

Z_cur = get_Z(shot, timebase_column);

% Next, get the programmed current centroid from the Plasma Control System
% (PCS) tree.  This is somewhat complicated, since the PCS architecture
% consists of multiple programming "segments".  Each segment covers a
% specific time range during the shot.  Each segment has a specified start
% time, but the segments are not necessarily in sequential order in time.
% The segments can be "ON" (i.e. active) or "OFF" (i.e. inactive).  Only
% Z_cur programming from the "ON" segments should be read in.  Each segment
% contains programming and PID gains for up to 16 different parameters.
% (These are called "wires" in PCS.)  In each active segment, we want to
% read in the Z_cur programming from either one of the two wires that are
% named "ZCUR" and contain non-zero PID gains.  If no such wires exist in
% an active segment, then the Z_cur programming is set to NaN for that
% segment's time range.  The C-Mod PCS has a clocking rate of 1 kHz, so we
% will assemble the full shot's Z_cur programming onto a 1 kHz timebase
% spanning t=-4 s to t=+12 s, with the programming from each segment's time
% range appropriately mapped onto it.

% The first step is to determine which segments are ON, and what their
% respective start times are.

[shotopened, status] = mdsopen('pcs',shot);
if (mod(status,2) == 1);

% Find the active segments

  states = mdsvalue('getnci("\\top.seg_*:seg_num", "state")');
  active_segments = find(states == 0);   % 0 --> "ON",  1 --> "OFF"
  nsegs = length(active_segments);

% Get the start time for each active segment

  start_times = single(zeros(nsegs, 1));
  for iseg = 1:nsegs;
    nodename = ['\top.seg_'  num2str(active_segments(iseg),'%02i')  ...
                ':start_time'];
    start_times(iseg) = single(mdsvalue(nodename));
  end;

% Now sort the start times into sequential order, and apply the same
% ordering to the list of active segments

  [start_times, ordered_indices] = sort(start_times);
  active_segments = active_segments(ordered_indices);
  stop_times = circshift(start_times, -1);
  stop_times(end) = 12.383;

% Create the 1 kHz PCS timebase.  Each segment's Ip programming will be
% mapped onto this timebase.

  pcstime = transpose(single([-4 : 0.001 : 12.383]));
  Z_prog = NaN(size(pcstime));
  Z_prog_temp = NaN(size(pcstime));

% For each active segment, find the first wire for Z_cur control, and
% determine if it has non-zero PID gains.  If it exists, interpolate its
% Z_cur programming onto the PCS timebase.  Finally, use only the
% interpolated Z_cur programming values that are within the specified start
% and stop times for each segment.

  for iseg = 1:nsegs;
    seg = active_segments(iseg);
    for iwire = 1:16;
      nodename = ['.seg_'  num2str(seg,'%02i')  ...
                  ':P_'  num2str(iwire,'%02i')  ':name'];
      wirename = mdsvalue(nodename);
      if (strcmpi(wirename, 'ZCUR'));
        nodename = ['.seg_'  num2str(seg,'%02i')  ...
                    ':P_'  num2str(iwire,'%02i')  ':pid_gains'];
        pid_gains = mdsvalue(nodename);
        if (~isempty(find(pid_gains > 0)));
          nodename = ['.seg_'  num2str(seg,'%02i')  ...
                      ':P_'  num2str(iwire,'%02i')];
          signal = mdsvalue(nodename);
          sigtime = mdsvalue(['dim_of(' nodename ')']);
          Z_prog_temp = interp1(sigtime, signal, pcstime, 'linear', ...
            signal(end));
          segment_time_indices = find(pcstime >= start_times(iseg) & ...
                                      pcstime <=  stop_times(iseg));
          Z_prog(segment_time_indices) = Z_prog_temp(segment_time_indices);
        end;
        break;
      end;
    end;
  end;
  mdsclose;
end;

% Interpolate Z_prog onto the requested timebase

Z_prog = interp1(pcstime, Z_prog, timebase_column, 'linear', Z_prog(end));

% Finally, calculate the Z error

Z_error = Z_prog - Z_cur;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Z_cur   = transpose(Z_cur);
  Z_prog  = transpose(Z_prog);
  Z_error = transpose(Z_error);
end;

end