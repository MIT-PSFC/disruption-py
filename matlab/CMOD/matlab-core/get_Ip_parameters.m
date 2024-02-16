function [ip, ip_prog, ip_error, dip_dt, dipprog_dt] = ...
  get_Ip_parameters(shot, timebase);

% This routine calculates Ip_error = (Ip - Ip_programmed), i.e. how much
% the actual plasma current differs from the requested current.  It
% linearly interpolates both the programmed and measured plasma currents
% onto the given timebase.  The routine also calculates the time
% derivatives of the measured Ip and programmed Ip.  The time derivatives
% are useful for discriminating between rampup, flattop, and rampdown.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ip = measured plasma current [A]
%   ip_prog = programmed (i.e. requested) plasma current [A]
%   ip_error = (Ip - Ip_programmed) [A]
%   dip_dt = d(Ip)/dt [A/s]
%   dipprog_dt = d(Ip_programmed)/dt [A/s]
%
% NOTE: Given ip, ip_prog, and ip_error, we can trivially calculate the
% relative errors, i.e. ip_error/ip and ip_error/ip_prog.  We may add one
% or both of these to the output parameter list at some point in the
% future.
%
% Authors: Alex Tinguely and Robert Granetz   Oct 2015
% Modified by:
% Ryan Sweeney, 18 Jan. 2019
% - added an else statement when checking the "status" of opening the mdstree.
		      
% 
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

ip         = NaN(length(timebase), 1);
ip_prog    = NaN(length(timebase), 1);
ip_error   = NaN(length(timebase), 1);
dip_dt     = NaN(length(timebase), 1);
dipprog_dt = NaN(length(timebase), 1);

% Next, read in the measured plasma current, Ip, from the magnetics tree.
% If Ip is successfully read in, then calculate dI/dt, and then interpolate
% both signals onto the requested timebase.

[shotopened, status] = mdsopen('magnetics',shot);
if (mod(status,2) == 1);
  ip = mdsvalue('\ip');
  [magtime, status] = mdsvalue('dim_of(\ip)');
  if (mod(status,2) == 1);
    dip_dt = gradient(ip, magtime);
    ip = interp1(magtime, ip, timebase_column, 'linear');
    dip_dt = interp1(magtime, dip_dt, timebase_column, 'linear');
  end;
  mdsclose;
end;

% Next, get the programmed plasma current from the Plasma Control System
% (PCS) tree.  This is somewhat complicated, since the PCS architecture
% consists of multiple programming "segments".  Each segment covers a
% specific time range during the shot.  Each segment has a specified start
% time, but the segments are not necessarily in sequential order in time.
% The segments can be "ON" (i.e. active) or "OFF" (i.e. inactive).  Only Ip
% programming from the "ON" segments should be read in.  Each segment
% contains programming and PID gains for up to 16 different parameters.
% (These are called "wires" in PCS.)  In each active segment, we want to
% read in the Ip programming from the wire that is named "IP" and contains
% non-zero PID gains.  If no such wire exists in an active segment, then
% the Ip programming is set to NaN for that segment's time range.  The
% C-Mod PCS has a clocking rate of 1 kHz, so we will assemble the full
% shot's Ip programming onto a 1 kHz timebase spanning t=-4 s to t=+12 s,
% with the programming from each segment's time range appropriately mapped
% onto it.

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
  ip_prog = NaN(size(pcstime));
  ip_prog_temp = NaN(size(pcstime));

% For each active segment, find the wire for Ip control, and determine if
% it has non-zero PID gains.  If it exists, interpolate its Ip programming
% onto the PCS timebase.  Finally, use only the interpolated Ip programming
% values that are within the specified start and stop times for each
% segment.

  for iseg = 1:nsegs;
    seg = active_segments(iseg);
    for iwire = 16 : -1 : 1;  % The Ip wire is usually #16, so start with that
      nodename = ['.seg_'  num2str(seg,'%02i')  ...
                  ':P_'  num2str(iwire,'%02i')  ':name'];
      wirename = mdsvalue(nodename);
      if (strcmpi(wirename, 'IP'));
        nodename = ['.seg_'  num2str(seg,'%02i')  ...
                    ':P_'  num2str(iwire,'%02i')  ':pid_gains'];
        pid_gains = mdsvalue(nodename);
        if (~isempty(find(pid_gains > 0)));
          nodename = ['.seg_'  num2str(seg,'%02i')  ...
                      ':P_'  num2str(iwire,'%02i')];
          signal = mdsvalue(nodename);
          sigtime = mdsvalue(['dim_of(' nodename ')']);
          ip_prog_temp = interp1(sigtime, signal, pcstime, 'linear', ...
            signal(end));
          segment_time_indices = find(pcstime >= start_times(iseg) & ...
                                      pcstime <=  stop_times(iseg));
          ip_prog(segment_time_indices) = ip_prog_temp(segment_time_indices);
        end;
        break;
      end;
    end;
  end;
  mdsclose;

% Calculate time derivative of ip_prog

  dipprog_dt = gradient(ip_prog, pcstime);

% RMS - added the else statement here to treat the case
% where opening the mdstree is unsuccessful. 			
else
  return
end;

% Interpolate ip_prog and dipprog_dt onto the requested timebase

ip_prog = interp1(pcstime, ip_prog, timebase_column, 'linear', ip_prog(end));
dipprog_dt = interp1(pcstime, dipprog_dt, timebase_column, 'linear', 0);

ip_error = ip - ip_prog; % calculate the ip_error

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  ip         = transpose(ip);
  ip_prog    = transpose(ip_prog);
  ip_error   = transpose(ip_error);
  dip_dt     = transpose(dip_dt);
  dipprog_dt = transpose(dipprog_dt);
end;

end
