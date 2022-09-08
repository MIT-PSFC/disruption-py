function [beta_N, beta_p, beta_p_dot, kappa, upper_gap, lower_gap, ...
  li, li_dot, q0, qstar, q95, V_loop_efit, Wmhd, dWmhd_dt, ...
  ssep, n_over_ncrit] = get_EFIT_parameters(shot, timebase);

% This function obtains the EFIT parameters, and also calculates the time
% derivatives of a few of them.
%
% Inputs:
%   shot = shot number
%   timebase = array of times at which to evaluate the parameters
%
% Outputs:
%   beta_N = normalized beta
%   beta_p = beta_poloidal
%   beta_p_dot = time derivative of beta_p [/s]
%   kappa = elongation
%   upper_gap = upper gap [m]
%   lower_gap = lower gap [m]
%   li = internal inductance
%   li_dot = time derivative of li [/s]
%   q0 = safety factor at center
%   qstar = q* = cylindrical safety factor
%   q95 = edge safety factor
%   V_loop_efit = loop voltage [V]
%   Wmhd = diamagnetic energy, [J]
%   Wmhd_dot = time derivative of Wmhd, [J/s]
%   ssep = distance on midplane between 1st and 2nd separatrices [m]
%   n_over_ncrit = vertical stability criterion (EFIT name: xnnc)
%
% Author: Alex Tinguely 2015-10-07 (Dan Brunner's 30th Birthday!)
% Some modifications by Robert Granetz 2015/10/20

% Initialize all outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

beta_N = NaN(size(timebase));
beta_p = NaN(size(timebase));
beta_p_dot = NaN(size(timebase));
elong = NaN(size(timebase));
upper_gap = NaN(size(timebase));
lower_gap = NaN(size(timebase));
li = NaN(size(timebase));
li_dot = NaN(size(timebase));
q0 = NaN(size(timebase));
qstar = NaN(size(timebase));
q95 = NaN(size(timebase));
V_loop_efit = NaN(size(timebase));
Wmhd = NaN(size(timebase));
ssep = NaN(size(timebase));
n_over_ncrit = NaN(size(timebase));

% All the signals used inside this routine are column vectors, since
% mdsvalue returns column vectors.  But the input parameter, 'timebase',
% may be either a column or row vector.  Therefore, in order to continue
% with our column-based calculations, we need to create a copy of
% 'timebase' that is guaranteed to be a column vector.  Just before
% returning from this routine, all output vectors will be converted to the
% same shape (i.e. column or row vector) as the 'timebase' input
% parameter.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% EFIT data on the standard timebase are in the ANALYSIS tree.  For shots
% that disrupted, these data are superceded by EFIT data in the EFIT18
% tree, which includes data at a high sampling rate prior to the disruption
% time, in addition to the standard times.  This routine does not know, a
% priori, which tree it should open.  Therefore, I will first try to open
% the EFIT18 tree, and if this fails, then I'll open the ANALYSIS tree.

[shotopened, status] = mdsopen('efit18', shot);
if (mod(status,2) == 0);
  [shotopened, status] = mdsopen('analysis', shot);
  if (mod(status,2)==0);
    return;
  end;
end;

% Read in EFIT timebase

[efittime, status] = mdsvalue('\efit_aeqdsk:time');
if (mod(status,2)==0 || length(efittime) == 1);
  return;
end;

beta_N = mdsvalue('\efit_aeqdsk:betan');
beta_p = mdsvalue('\efit_aeqdsk:betap');
kappa = mdsvalue('\efit_aeqdsk:eout');
li = mdsvalue('\efit_aeqdsk:li');
upper_gap = mdsvalue('\efit_aeqdsk:otop')/100; % meters
lower_gap = mdsvalue('\efit_aeqdsk:obott')/100; % meters
q0 = mdsvalue('\efit_aeqdsk:q0');
qstar = mdsvalue('\efit_aeqdsk:qstar');
q95 = mdsvalue('\efit_aeqdsk:q95');
V_loop_efit = mdsvalue('\efit_aeqdsk:vloopt');
Wmhd = mdsvalue('\efit_aeqdsk:wplasm');
ssep = mdsvalue('\efit_aeqdsk:ssep')/100; % meters
n_over_ncrit =mdsvalue('-\efit_aeqdsk:xnnc');

mdsclose;

beta_p_dot = gradient(beta_p, efittime);
li_dot = gradient(li, efittime);
dWmhd_dt = gradient(Wmhd, efittime);

% Interpolate data onto the requested timebase

beta_N = interp1(efittime, beta_N, timebase_column, 'linear');
beta_p = interp1(efittime, beta_p, timebase_column, 'linear');
kappa = interp1(efittime, kappa, timebase_column, 'linear');
li = interp1(efittime, li, timebase_column, 'linear');
upper_gap = interp1(efittime, upper_gap, timebase_column, 'linear');
lower_gap = interp1(efittime, lower_gap, timebase_column, 'linear');
q0= interp1(efittime, q0, timebase_column, 'linear');
qstar = interp1(efittime, qstar, timebase_column, 'linear');
q95 = interp1(efittime, q95, timebase_column, 'linear');
V_loop_efit = interp1(efittime, V_loop_efit, timebase_column, 'linear');
Wmhd = interp1(efittime, Wmhd, timebase_column, 'linear');
beta_p_dot = interp1(efittime, beta_p_dot, timebase_column, 'linear');
li_dot = interp1(efittime, li_dot, timebase_column, 'linear');
dWmhd_dt = interp1(efittime, dWmhd_dt, timebase_column, 'linear');
ssep = interp1(efittime, ssep, timebase_column, 'linear');
n_over_ncrit = interp1(efittime, n_over_ncrit, timebase_column, 'linear');

if (size(timebase,2) > 1);    % If "timebase" is a row vector, then convert
  beta_N = transpose(beta_N); % all the outputs back to row vectors
  beta_p = transpose(beta_p);
  kappa = transpose(kappa);
  li = transpose(li);
  upper_gap = transpose(upper_gap);
  lower_gap = transpose(lower_gap);
  q0 = transpose(q0);
  qstar = transpose(qstar);
  q95 = transpose(q95);
  V_loop_efit = transpose(V_loop_efit);
  Wmhd = transpose(Wmhd);
  beta_p_dot = transpose(beta_p_dot);
  li_dot = transpose(li_dot);
  dWmhd_dt = transpose(dWmhd_dt);
  ssep = transpose(ssep);
  n_over_ncrit = transpose(n_over_ncrit);
end;

end