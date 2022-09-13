function [Z_error, Z_prog, Z_cur, v_z, z_times_v_z] = ...
    get_Z_parameters(shot, timebase);

% This script was adapted from get_Z_error.m.
%
% The purpose of this script is to read in the values of Z_error and Z_prog
% from the plasma control system (PCS). Z_prog is the programmed vertical
% position of the plasma current centroid, and Z_error is the difference
% between the actual position and that requested (Z_error = Z_cur -
% Z_prog). Thus, the actual (estimated) position, Z_cur, can be calculated.
% And the vertical velocity, v_z, can be taken from the time derivative,
% and the product z_times_v_z ( = Z_cur * v_z) is also calculated. All of
% these values are then linearly interpolated over the given timebase.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values [s]
%
% Outputs:
%   Z_error from PCS [m]
%   Z_prog = programmed (i.e. requested) vertical position of the current
%     centroid [m]
%   Z_cur = actual vertical position of the current centroid, calculated
%         = Z_cur - Z_prog [m]
%   v_z = time derivative of Z_cur [m/s]
%   z_times_v_z = Z_cur * v_z [m^2/s]
%
% Original Author: Alex Tinguely 2015-09-09
% Rewritten by Robert Granetz 2016/01/29
% Updated by Alex Tinguely 2016-05-02
% Updated by Robert Granetz 2017/12/22 to use Ip from magnetics tree instead
%   of the PCS tree for shots prior to 2015.  This is because a factor
%   required to convert the PCS Ip signal to units of amperes is apparently
%   not available in the PCS tree for shots prior to the 2015 run campaign
% Updated by Robert Granetz 2020/06/26 to invert the signs of Z_cur and
%   Z_error, so they agree with my intended convention.  Apparently the PCS
%   signals are inverted from what I expected.  This new sign change occurs
%   at the very of this program.  I did not change anything related to
%   Z_prog, v_z, or z_times_v_z.

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
v_z     = NaN(length(timebase),1);
z_times_v_z   = NaN(length(timebase),1);

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
% segment's time range. The C-Mod PCS has a clocking rate of 1 kHz, so we
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

% Create a variable wire_index to be used when looking for Z_error later.

  wire_index = 0;
  
% For each active segment, find the first wire for Z_cur control, and
% determine if it has non-zero PID gains.  If it exists, interpolate its
% Z_cur programming onto the DPCS timebase.  Finally, use only the
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
          wire_index = iwire; % save this for Z_error below  
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

% Now we call Z_error from the Hybrid tree.

[shotopened, status] = mdsopen('hybrid',shot);
if (mod(status,2) == 1);
    
% Read in A_OUT, which is a 16xN matrix of the errors for *all* 16 wires for
% *all* of the segments. Note that DPCS time is usually taken at 10kHz.

    error_matrix = mdsvalue('\top.hardware.dpcs.signals:a_out');
    dpcstime = mdsvalue('dim_of(\top.hardware.dpcs.signals:a_out,1)'); % s

% The value of Z_error we read is not in the units we want. It must be
% *divided* by a factor AND *divided* by the plasma current.

    Z_error_without_factor_and_ip = error_matrix(wire_index,:);
    Z_error_without_ip = NaN(size(Z_error_without_factor_and_ip));
    Z_error = NaN(size(Z_error_without_factor_and_ip));

% Also, it turns out that different segments have different factors. So we
% search through the active segments (determined above), find the factors,
% and *divide* by the factor only for the times in the active segment (as
% determined from start_times and stop_times.
    
    for isegs = 1:nsegs % for each active segment
        seg = active_segments(iseg);
        nodename = ['\dpcs::top.seg_' num2str(seg,'%02i') ...
            ':p_' num2str(wire_index,'%02i') ':predictor:factor'];
        Z_factor = mdsvalue(nodename);
        segment_time_indices = find(dpcstime >= start_times(iseg) & ...
                                      dpcstime <=  stop_times(iseg));
        Z_error_without_ip(segment_time_indices) = ...
            Z_error_without_factor_and_ip(segment_time_indices)/Z_factor; % [A*m]  
    end
    
% Next we grab ip, which comes from a_in:input_056. This also requires 
% *multiplication* by a factor.

% NOTE that I can't get the following ip_without_factor to work for shots
% before 2015.

    if (shot >= 1150101000);
      ip_without_factor = transpose(...
          mdsvalue('\hybrid::top.hardware.dpcs.signals.a_in:input_056'));
%     iptime = mdsvalue(['dim_of(\hybrid::top.hardware.dpcs.signals' ...
%         '.a_in:input_056)']);
      
      ip_factor = mdsvalue(['\hybrid::top.dpcs_config.inputs:'...
          'input_056:p_to_v_expr']);
      
      ip = ip_without_factor*ip_factor; % [A]
    else;
     mdsopen('magnetics', shot);
     ip = mdsvalue('\ip');
     iptime = mdsvalue('dim_of(\ip)');
     mdsclose;
     ip = transpose(interp1(iptime, ip, dpcstime, 'linear'));
    end;
      
    Z_error = Z_error_without_ip./ip; % [m] 
    
end

% We calculate Z_cur = Z_prog + Z_error. Currently, they could be on
% different timebases, but ultimately the original Z_prog was just a few
% points, so we interpolate Z_prog on the faster Z_error timebase.

Z_prog_dpcs = transpose(interp1(pcstime, Z_prog, dpcstime));
Z_cur = Z_prog_dpcs + Z_error; % m

% We calculate the vertical velocity v_z

v_z = gradient(Z_cur, dpcstime); % m/s
z_times_v_z = Z_cur.*v_z; % m^2/s

% Interpolate all signals onto the requested timebase

Z_prog = interp1(pcstime, Z_prog, timebase_column, 'linear', Z_prog(end));
Z_error = interp1(dpcstime, Z_error, timebase_column, 'linear',...
    Z_error(end));
Z_cur = interp1(dpcstime, Z_cur, timebase_column, 'linear', Z_cur(end));
v_z = interp1(dpcstime, v_z, timebase_column, 'linear', v_z(end));
z_times_v_z = interp1(dpcstime, z_times_v_z, timebase_column, 'linear',...
    z_times_v_z(end));

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Z_cur   = transpose(Z_cur);
  Z_prog  = transpose(Z_prog);
  Z_error = transpose(Z_error);
  v_z     = transpose(v_z);
  z_times_v_z = transpose(z_times_v_z);
end;

Z_cur = -Z_cur;
Z_error = -Z_error;

end