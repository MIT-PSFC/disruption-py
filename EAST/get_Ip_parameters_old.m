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
% Author: Robert Granetz   Dec 2015

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

ip         = NaN(length(timebase), 1);
ip_prog    = NaN(length(timebase), 1);
ip_error   = NaN(length(timebase), 1);
dip_dt     = NaN(length(timebase), 1);
dipprog_dt = NaN(length(timebase), 1);

% Next, read in the measured plasma current, Ip.  There are several
% different measurements of Ip: IPE, IPG, IPM (all in the EAST tree), and
% PCRL01 (in the PCS_EAST tree).  At various times in the history of EAST,
% there are been problems with all of these measurements, such as broken
% sensors, inverted signals, and shifted timebases.  I think the most
% reliable one is PCRL01, which is the one used by the Plasma Control
% System (PCS) for feedback control.  So that is the one I will use for the
% disruption warning database.

% Get \PCRL01 plasma current signal from the PCS_EAST tree.  If there is no
% PCS_EAST tree, or no PCRL01 data, then this is not a valid plasma shot
% (as per guidance from Qian Jinping).

[shotopened, status] = mdsopen('pcs_east', double(shot));
if (mod(status,2) == 1);
  [ip, status] = mdsvalue('\pcrl01');  % Read in Ip data (amps)
  if (mod(status,2) == 1);             % If successful, continue
    ip_time = mdsvalue('dim_of(\pcrl01)');

% For shots before year 2014, the PCRL01 timebase needs to be shifted by
% 17.0 ms 

    if (shot < 44432);
      ip_time = ip_time - 0.0170;
    end;

% High-frequency noise spikes on some shots can cause a problem with the
% time derivative and other computations.  Use a median filter to reduce
% the problem.

    ip = medfilt1(ip, 5); % Remove noise spikes with median filter

% Subtract baseline offset

    baseindices=find(ip_time <= -5.8); % time before any PF supplies turn on
    if (numel(baseindices) > 0);
      baseline=sum(ip(baseindices))/length(baseindices);
      ip = ip - baseline;
    end;

% Now calculate dIp/dt, and then interpolate both signals onto the
% requested timebase.

    dip_dt = gradient(ip, ip_time);
    ip = interp1(ip_time, ip, timebase_column, 'linear');
    dip_dt = interp1(ip_time, dip_dt, timebase_column, 'linear');
  end;

% Next, get the programmed plasma current.

  [ip_prog, status] = mdsvalue('\lmtipref');
  if (mod(status,2) == 1);             % If successful, continue
    ip_prog = ip_prog * 1.e6; % convert from MA to A
    ip_prog_time = mdsvalue('dim_of(\lmtipref)');

% For shots before year 2014, the LMTIPREF timebase needs to be shifted by
% 17.0 ms 

    if (shot < 44432);
      ip_prog_time = ip_prog_time - 0.0170;
    end;

% Calculate time derivative of ip_prog

    dipprog_dt = gradient(ip_prog, ip_prog_time);

% Interpolate ip_prog and dipprog_dt onto the requested timebase

    ip_prog = interp1(ip_prog_time, ip_prog, timebase_column, 'linear');
    dipprog_dt = interp1(ip_prog_time, dipprog_dt, timebase_column,'linear');

    ip_error = ip - ip_prog; % calculate the ip_error
  end;
  mdsclose;
end;

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
