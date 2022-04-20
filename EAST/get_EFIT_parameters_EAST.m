function [beta_n, beta_p, dbetap_dt, kappa, upper_gap, lower_gap, ...
  li, dli_dt, q0, qstar, q95, Wmhd, dWmhd_dt] = ...
  get_EFIT_parameters_EAST(shot, timebase);

% This function obtains the EFIT parameters, and also calculates the time
% derivatives of a few of them.
%
% Inputs:
%   shot = shot number
%   timebase = array of times at which to evaluate the parameters
%
% Outputs:
%   beta_n = normalized beta
%   beta_p = beta poloidal
%   dbetap_dt = time derivative of beta_p [/s]
%   kappa = elongation
%   upper_gap = upper gap [m]
%   lower_gap = lower gap [m]
%   li = internal inductance
%   dli_dt = time derivative of li [/s]
%   q0 = safety factor at center
%   qstar = q* = cylindrical safety factor
%   q95 = edge safety factor
%   Wmhd = diamagnetic energy, [J]
%   dWmhd_dt = time derivative of Wmhd, [J/s]
%  
%   Note that:
%   (1) No V_loop data in "efit_east" tree on EAST
%   (2) Some EFIT node names are different from C-Mod
%       'kappa'  instead of 'eout'
%       'gaptop' instead of 'otop'
%       'gapbot' instead of 'obott'  
%
% Author: Modified on EAST by Wang Bo 2015/12/10
%         Alex Tinguely 2015-10-07 
%         Some modifications on EAST by Robert Granetz 2016/04/05
%         2017/04/10 - R. Granetz; modified so that program tries to the
%           standard efit tree, 'efit_east', in the case when our special
%           efit tree, 'efit18', does not exist.
%
% Revision/history:
%  2017/12/13 - RSG; Use chisq to determine invalid time slices
%  2019/12/22 - RSG & Jinxiang Zhu; add in routine to explicitly calculate
%               upper and lower gaps from EFIT information on plasma
%               boundary and first wall geometry.  We needed to do this
%               because the EAST version of EFIT does not output the upper
%               and lower gaps directly.

% Initialize all outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

beta_n     = NaN(size(timebase));
beta_p     = NaN(size(timebase));
dbetap_dt  = NaN(size(timebase));
kappa      = NaN(size(timebase));
upper_gap  = NaN(size(timebase));
lower_gap  = NaN(size(timebase));
li         = NaN(size(timebase));
dli_dt     = NaN(size(timebase));
q0         = NaN(size(timebase));
qstar      = NaN(size(timebase));
q95        = NaN(size(timebase));
Wmhd       = NaN(size(timebase));
dWmhd_dt   = NaN(size(timebase));

% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, in Matlab the "mdsvalue" routine returns
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

% Open EFIT18 tree.  If error, try opening the standard EFIT tree,
% 'EFIT_EAST'.  If that also fails, then exit.

[shotopened, status] = mdsopen('efit18', double(shot));
if (mod(status,2) == 0);
  [shotopened, status] = mdsopen('efit_east', double(shot));
  if (mod(status,2) == 0);
    return;
  end;
end;
% Read in the EFIT timebase.  If error, exit routine.

[efittime, status] = mdsvalue('\efit_aeqdsk:atime');
[efittime, unique_indices, ~] = unique(efittime);
if (mod(status,2)==0 || length(efittime) <= 1);
  return;
end;

% Read in the EFIT signals.  Assume that if the first signal is okay, all
% the signals will all be okay.  And if the first signal has an error, they
% will all have errors, and we should exit the routine.

[beta_n, status] = mdsvalue('\efit_aeqdsk:betan');
if (mod(status,2)==0);
  beta_n = NaN(size(timebase));
  return;
end;
beta_p = mdsvalue('\efit_aeqdsk:betap'); 
kappa = mdsvalue('\efit_aeqdsk:kappa'); % different from C-Mod
li = mdsvalue('\efit_aeqdsk:li');
q0 = mdsvalue('\efit_aeqdsk:q0');
qstar = mdsvalue('\efit_aeqdsk:qstar');
q95 = mdsvalue('\efit_aeqdsk:q95');
Wmhd = mdsvalue('\efit_aeqdsk:wmhd');
chisq = mdsvalue('\efit_aeqdsk:chisq'); % Use chisq to determine which time
                                        % slices are invalid.
%upper_gap= mdsvalue('\efit_aeqdsk:gaptop')/100; % meters, different from C-Mod
%lower_gap= mdsvalue('\efit_aeqdsk:gapbot')/100; % meters, different from C-Mod
mdsclose;

[upper_gap, lower_gap] = get_EFIT_gaps_EAST(double(shot)); % Jinxiang Zhu's
                                                           % routine
% Note: Jinxiang's routine does an mdsopen and mdsclose of the EFIT tree,
% so it needs to be called after all the other EFIT information has been
% obtained.

beta_n = beta_n(unique_indices);
beta_p = beta_p(unique_indices);
kappa = kappa(unique_indices);
li = li(unique_indices);
upper_gap = upper_gap(unique_indices);
lower_gap = lower_gap(unique_indices);
q0 = q0(unique_indices);
qstar = qstar(unique_indices);
q95 = q95(unique_indices);
Wmhd = Wmhd(unique_indices);
chisq = chisq(unique_indices);

% EFIT reconstructions are sometimes invalid, particularly when very close
% to a disruption.  There are a number of EFIT parameters that can indicate
% invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
% 'chisq' to determine which time slices should be excluded from our
% disruption warning database.

invalid_indices = find(chisq > 50);

beta_n(invalid_indices) = NaN;
beta_p(invalid_indices) = NaN;
kappa(invalid_indices) = NaN;
li(invalid_indices) = NaN;
upper_gap(invalid_indices) = NaN;
lower_gap(invalid_indices) = NaN;
q0(invalid_indices) = NaN;
qstar(invalid_indices) = NaN;
q95(invalid_indices) = NaN;
Wmhd(invalid_indices) = NaN;

% Take time derivative of some of the EFIT signals

dbetap_dt = gradient(beta_p, efittime);
dli_dt = gradient(li, efittime);
dWmhd_dt = gradient(Wmhd, efittime);

% Interpolate data onto the requested timebase.  This routine used to have
% 'interp1' extrapolate, using the value 'zero' for the extrapolation.  We
% should not do that for EFIT parameters.  If extrapolation is needed, the
% extrapolated values should be 'NaN' for EFIT parameters.

beta_n    = interp1(efittime, beta_n,    timebase_column);
beta_p    = interp1(efittime, beta_p,    timebase_column);
kappa     = interp1(efittime, kappa,     timebase_column);
li        = interp1(efittime, li,        timebase_column);
upper_gap = interp1(efittime, upper_gap, timebase_column);
lower_gap = interp1(efittime, lower_gap, timebase_column);
q0        = interp1(efittime, q0,        timebase_column);
qstar     = interp1(efittime, qstar,     timebase_column);
q95       = interp1(efittime, q95,       timebase_column);
Wmhd      = interp1(efittime, Wmhd,      timebase_column);
dbetap_dt = interp1(efittime, dbetap_dt, timebase_column);
dli_dt    = interp1(efittime, dli_dt,    timebase_column);
dWmhd_dt  = interp1(efittime, dWmhd_dt,  timebase_column);

% If "timebase" is a row vector, then convert all the outputs back to row
% vectors before exiting this routine.

if (size(timebase,2) > 1);
  beta_n    = transpose(beta_n);
  beta_p    = transpose(beta_p);
  kappa     = transpose(kappa);
  li        = transpose(li);
  upper_gap = transpose(upper_gap);
  lower_gap = transpose(lower_gap);
  q0        = transpose(q0);
  qstar     = transpose(qstar);
  q95       = transpose(q95);
  Wmhd      = transpose(Wmhd);
  dbetap_dt = transpose(dbetap_dt);
  dli_dt    = transpose(dli_dt);
  dWmhd_dt  = transpose(dWmhd_dt);
end;

end
