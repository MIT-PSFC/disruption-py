function   [ip_RT, ip_prog_RT, ip_error_RT, dip_dt_RT, dipprog_dt_RT, ...
    power_supply_railed] = get_Ip_parameters_RT(shot, timebase);

% This routine calculates Ip_error = (Ip - Ip_programmed), i.e. how much
% the actual plasma current differs from the requested current, based on
% the real-time signals available from PCS.  It linearly interpolates both
% the programmed and measured plasma currents onto the given timebase.  The
% routine also calculates the time derivatives of the measured Ip and
% programmed Ip.  The time derivatives are useful for discriminating
% between rampup, flattop, and rampdown.  Finally, it also looks at the
% signal 'epsoff' to determine if one of the E-coil power supplies has
% railed.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ip_RT = measured plasma current [A]
%   ip_prog_RT = programmed (i.e. requested) plasma current [A]
%   ip_error_RT = (Ip - Ip_programmed), unless Ip feedback is disabled [A]
%   dip_dt_RT = d(Ip)/dt [A/s]
%   dipprog_dt_RT = d(Ip_programmed)/dt [A/s]
%
% NOTE: Given ip, ip_prog, and ip_error, we can trivially calculate the
% relative errors, i.e. ip_error/ip and ip_error/ip_prog.  We may add one
% or both of these to the output parameter list at some point in the
% future.
%
% Authors: Robert Granetz   Aug 2017
%   Modified from get_Ip_parameters.m

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

ip_RT         = NaN(length(timebase), 1);
ip_prog_RT    = NaN(length(timebase), 1);
ip_error_RT   = NaN(length(timebase), 1);
dip_dt_RT     = NaN(length(timebase), 1);
dipprog_dt_RT = NaN(length(timebase), 1);
power_supply_railed = NaN(length(timebase), 1);

% Read in the measured plasma current using the real-time signal available
% in PCS.  If Ip is successfully read in, then calculate dI/dt, and then
% interpolate both signals onto the requested timebase.

[iptime, status] = mdsvalue(['dim_of(ptdata("ipsip", ' num2str(shot) '))']);
if (mod(status,2) == 1);
  ip_RT = mdsvalue(['ptdata("ipsip", ' num2str(shot) ')']);
else;
  [iptime, status] = mdsvalue(['dim_of(ptdata("ipspr15V", ' ...
    num2str(shot) '))']);
  if (mod(status,2) == 1);
    ip_RT = mdsvalue(['ptdata("ipspr15V", ' num2str(shot) ')']) / 2;
  end;
end;

if (mod(status,2) == 1);
  ip_RT = ip_RT * 1.e6; % convert to amperes
  iptime = iptime/1.e3; % convert ms to s
  dip_dt_RT = gradient(ip_RT, iptime);
  ip_RT = interp1(iptime, ip_RT, timebase_column, 'linear');
  dip_dt_RT = interp1(iptime, dip_dt_RT, timebase_column, 'linear');
end;

% Next, get the programmed plasma current.

[ip_prog_time, status] = mdsvalue(['dim_of(ptdata("ipsiptargt", ' ...
                                                   num2str(shot) '))']);
if (mod(status,2) == 1);
  ip_prog_RT = mdsvalue(['ptdata("ipsiptargt", ' num2str(shot) ')']) * 0.5;
  ip_prog_RT = ip_prog_RT * 1.e6; % convert to amperers
  ip_prog_time = ip_prog_time/1.e3; % convert ms to s
  polarity = unique(mdsvalue(['ptdata("iptdirect", ' ...
                                                   num2str(shot) ')']));
  if (length(polarity) ~= 1);
    fprintf(1,'ERROR: polarity of Ip target is not constant\n');
    polarity = polarity(1);
  end;
  ip_prog_RT = ip_prog_RT * polarity;
  dipprog_dt_RT = ...
    gradient(ip_prog_RT, ip_prog_time); % Calculate time derivative
% Interpolate ip_prog_RT and dipprog_dt_RT onto the requested timebase
  ip_prog_RT = interp1(ip_prog_time, ip_prog_RT, timebase_column, 'linear');
  dipprog_dt_RT = ...
    interp1(ip_prog_time, dipprog_dt_RT, timebase_column, 'linear');
end;

% Now get 'ip_error_RT' directly from PCS

[ip_error_time, status] = mdsvalue(['dim_of(ptdata("ipeecoil", ' ...
                                                   num2str(shot) '))']);
if (mod(status,2) == 1);
  ip_error_time = ip_error_time/1.e3; % convert to seconds
  ip_error_RT = mdsvalue(['ptdata("ipeecoil", ' num2str(shot) ')']) * 0.5;
  ip_error_RT = ip_error_RT * 1.e6; % convert to amperes
  ip_error_RT = interp1(ip_error_time, ip_error_RT, timebase_column, 'linear');
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

indices_feedback_off = find(ipimode ~= 0 & ipimode ~= 3);
ip_error_RT(indices_feedback_off) = NaN;

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
end;

% Times at which power_supply_railed ~=0 (i.e. epsoff ~=0) mean that
% PCS feedback control of Ip is not being applied.  Therefore the
% 'ip_error' parameter is undefined for these times.

ps_railed_indices = find(power_supply_railed ~= 0);
ip_error_RT(ps_railed_indices) = NaN;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  ip_RT         = transpose(ip_RT);
  ip_prog_RT    = transpose(ip_prog_RT);
  ip_error_RT   = transpose(ip_error_RT);
  dip_dt_RT     = transpose(dip_dt_RT);
  dipprog_dt_RT = transpose(dipprog_dt_RT);
  power_supply_railed = transpose(power_supply_railed);
end;

end
