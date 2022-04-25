function P_ohm = get_P_ohm_EAST(shot, timebase);
%
% This function calculates the ohmic power, P_ohm.  The output vector,
% P_ohm, has the same shape as the input vector, "timebase", i.e. if
% "timebase" is a row vector, P_ohm will be a row vector; if "timebase" is
% a column vector, P_ohm will be a column vector.  Note that Matlab calls
% to mdsvalue return column vectors, so this routine works internally with
% column vectors.
%
% We use the following expression to calculate P_ohm:
%   P_ohm = Ip * V_resistive
%      where: V_resistive = V_loop - V_inductive
%                         = V_loop - L * dIp/dt
%      and L = L_internal = mu0 * R0 * li/2
%
%   pcs_east node '\pcvloop'        is used for V_loop,
%   efit18   node '\efit_aeqdsk:li' is used for li,
%   pcs_east node '\pcrl01'         is used for Ip
%
% and the time derivative of Ip is calculated using Matlab's "gradient"
% function with some boxcar smoothing.
%
% If the EFIT data or the magnetics Ip data is not available, then P_ohm is
% returned as an array of NaN's (Not-a-Number).
%
% Inputs:
%   timebase = times to calculate P_ohm at 
% Outputs:
%   P_ohm = ohmic power [W]
%
% Author: Modified on EAST by Wang Bo 2015/12/10
% Alex Tinguely 15-09-16
% Some modifications by Robert Granetz 2015/10/10
% More modifications to get this to work. by R. Granetz 2016/06/23

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

% Initialize

P_ohm = NaN(length(timebase), 1);

% Get loop voltage

[shotopened, status] = mdsopen('pcs_east', double(shot));
if (mod(status,2) == 1);
  [vloop, status] = mdsvalue('\pcvloop'); % volts
  if (mod(status,2) == 1);
    vloop_time = mdsvalue('dim_of(\pcvloop)'); % seconds
    vloop = interp1(vloop_time, vloop, timebase_column, 'linear', 0);
  else;
    vloop = NaN(length(timebase), 1);
  end;
else;
  vloop = NaN(length(timebase), 1);
end;

% Get internal inductance

[shotopened, status] = mdsopen('efit18',double(shot));
if (mod(status,2) == 1);
  [li_time, status] = mdsvalue('dim_of(\efit_aeqdsk:li)');
  if (mod(status,2) == 1 && length(li_time) > 1);
    li = mdsvalue('\efit_aeqdsk:li');
    [li_time, indx, ~] = unique(li_time);
    li = li(indx);
    li = interp1(li_time, li, timebase_column, 'linear', 0);
  else;
    li = ones(length(timebase), 1); % if error, use 1 instead of NaN
  end;
  mdsclose;
else;
  li = ones(length(timebase), 1); % if error, use 1 instead of NaN
end;
inductance = 4*pi*1.e-7 * 1.85 * li/2; % for EAST use R0 = 1.85 m
                                   
% Get plasma current from PCS

[shotopened, status] = mdsopen('pcs_east', double(shot));
if (mod(status,2) == 1);
  [ip, status] = mdsvalue('\pcrl01'); % amperes
  if (mod(status,2) == 1);
    ip_time = mdsvalue('dim_of(\pcrl01)');
    sign_ip = sign(sum(ip));
% calculate di/dt
    dipdt = gradient(ip, ip_time);
    dipdt_smoothed = smooth(dipdt, 11); % use 11-point boxcar smoothing
    ip = interp1(ip_time, ip, timebase_column, 'linear', 0);
    dipdt_smoothed = interp1(ip_time, dipdt_smoothed, ...
      timebase_column, 'linear', 0);
  else;
    dipdt_smoothed = NaN(length(timebase),1);
  end;
  mdsclose;
else;
  ip = NaN(length(timebase), 1);
  dipdt_smoothed = NaN(length(timebase), 1);
end;
 
v_inductive = -inductance .* dipdt_smoothed;
v_resistive = vloop - v_inductive;
v_resistive(find(v_resistive * sign_ip < 0)) = 0;
P_ohm = ip .* v_resistive;

% All of our calculations have used and produced column vectors, since
% mdsvalue returns column vectors.  If the input parameter, 'timebase', is
% a row vector, then make P_ohm a row vector.

if (size(timebase,2) > 1);
  P_ohm = transpose(P_ohm);
end;

end
