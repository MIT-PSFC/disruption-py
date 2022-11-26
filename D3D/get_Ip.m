function [ip, ip_prog, ip_error, dip_dt, dipprog_dt, power_supply_railed] = ...
  get_Ip_parameters(shot, timebase);

% This routine calculates Ip_error = (Ip - Ip_programmed), i.e. how much
% the actual plasma current differs from the requested current.  It
% linearly interpolates both the programmed and measured plasma currents
% onto the given timebase.  The routine also calculates the time
% derivatives of the measured Ip and programmed Ip.  The time derivatives
% are useful for discriminating between rampup, flattop, and rampdown.
% Finally, it also looks at the signal 'epsoff' to determine if one of the
% E-coil power supplies has railed.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ip = measured plasma current [A]
%   ip_prog = programmed (i.e. requested) plasma current [A]
%   ip_error = (Ip - Ip_programmed), unless Ip feedback is disabled [A]
%   dip_dt = d(Ip)/dt [A/s]
%   dipprog_dt = d(Ip_programmed)/dt [A/s]
%
% NOTE: Given ip, ip_prog, and ip_error, we can trivially calculate the
% relative errors, i.e. ip_error/ip and ip_error/ip_prog.  We may add one
% or both of these to the output parameter list at some point in the
% future.
%
% Authors: Alex Tinguely and Robert Granetz   Oct 2015

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
power_supply_railed = NaN(length(timebase), 1);

mdsconnect('atlas.gat.com');

[shotopened, status]=mdsopen('d3d', shot);
if (mod(status,2)==0);
% fprintf(1,'  Unable to open D3D tree for shot%7i\n', shot);
  if (size(timebase,2) > 1);
    ip         = transpose(ip);
    ip_prog    = transpose(ip_prog);
    ip_error   = transpose(ip_error);
    dip_dt     = transpose(dip_dt);
    dipprog_dt = transpose(dipprog_dt);
    power_supply_railed = transpose(power_supply_railed);
  end;
  return;
end;

% Next, read in the measured plasma current, Ip.
% If Ip is successfully read in, then calculate dI/dt, and then interpolate
% both signals onto the requested timebase.

ip = mdsvalue(['ptdata("ip", ' num2str(shot) ')']);
[iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' num2str(shot) '))']);
if (mod(status,2) == 1);
  iptime = iptime/1.e3; % convert ms to s
  dip_dt = gradient(ip, iptime);
  ip = interp1(iptime, ip, timebase_column, 'linear');
  dip_dt = interp1(iptime, dip_dt, timebase_column, 'linear');
else;
  ip = NaN(length(timebase), 1);
  dip_dt = NaN(length(timebase), 1);
end;

% Next, get the programmed plasma current.

ip_prog = mdsvalue(['ptdata("iptipp", ' num2str(shot) ')']);
[ip_prog_time, status] = mdsvalue(['dim_of(ptdata("iptipp", ' ...
                                                   num2str(shot) '))']);
if (mod(status,2) == 1);
  ip_prog_time = ip_prog_time/1.e3; % convert ms to s
  polarity = unique(mdsvalue(['ptdata("iptdirect", ' ...
                                                   num2str(shot) ')']));
  if (length(polarity) ~= 1);
    fprintf(1,'ERROR: polarity of Ip target is not constant\n');
    polarity = polarity(1);
  end;
  ip_prog = ip_prog * polarity;
  dipprog_dt = gradient(ip_prog, ip_prog_time); % Calculate time derivative
%   Interpolate ip_prog and dipprog_dt onto the requested timebase
  ip_prog = interp1(ip_prog_time, ip_prog, timebase_column, 'linear');
  dipprog_dt = interp1(ip_prog_time, dipprog_dt, timebase_column, 'linear');
else;
  ip_prog = NaN(length(timebase), 1);
  dipprog_dt = NaN(length(timebase), 1);
end;

% Now get the signal pointname 'ipimode'.  This PCS signal denotes whether
% or not PCS is actually feedback controlling the plasma current.  There
% are times when feedback of Ip is purposely turned off, such as during
% electron cyclotron current drive experiments.  Here is how to interpret
% the value of 'ipimode':
%  0: normal Ip feedback to E-coils supplies
%  3: almost normal Ip feedback, except that abs(Ip) > 2.5 MA
%  Anything else: not in normal Ip feedback mode.  In this case, the
% 'ip_prog' signal is irrelevant, and therefore 'ip_error' is not defined.

ipimode = mdsvalue(['ptdata("ipimode", ' num2str(shot) ')']);
[ipimode_time, status] = mdsvalue(['dim_of(ptdata("ipimode", ' ...
                                                   num2str(shot) '))']);
if (mod(status,2) == 1);
  ipimode_time = ipimode_time/1.e3;
  ipimode = interp1(ipimode_time, ipimode, timebase_column, 'linear');
else;
  ipimode = NaN(length(timebase), 1);
end;

% Now calculate 'ip_error':

ip_error = NaN(length(timebase), 1);
indices_feedback_on = find(ipimode == 0 | ipimode == 3);
ip_error(indices_feedback_on) = ip(indices_feedback_on) - ...
  ip_prog(indices_feedback_on);

% Finally, get 'epsoff' to determine if/when the E-coil power supplies have
% railed.

epsoff = mdsvalue(['ptdata("epsoff", ' num2str(shot) ')']);
[epsoff_time, status] = mdsvalue(['dim_of(ptdata("epsoff", ' ...
                                                  num2str(shot) '))']);
if (mod(status,2) == 1);
  epsoff_time = epsoff_time/1.e3;
  epsoff_time = epsoff_time + 0.001; % Avoid problem with simultaneity of
                                     % epsoff being triggered exactly on
				     % the last time sample
  epsoff = interp1(epsoff_time, epsoff, timebase_column, 'linear');
  power_supply_railed = zeros(length(timebase), 1);
  railed_indices = find(abs(epsoff) > 0.5);
  power_supply_railed(railed_indices) = 1;
else;
  power_supply_railed = NaN(length(timebase), 1);
end;

% Times at which power_supply_railed ~=0 (i.e. epsoff ~=0) mean that
% PCS feedback control of Ip is not being applied.  Therefore the
% 'ip_error' parameter is undefined for these times.

ps_railed_indices = find(power_supply_railed ~= 0);
ip_error(ps_railed_indices) = NaN;

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
  power_supply_railed = transpose(power_supply_railed);
end;

end
