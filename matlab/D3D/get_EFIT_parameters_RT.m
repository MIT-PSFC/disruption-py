function [beta_N_RT, beta_p_RT, dbeta_p_dt_RT, kappa_RT, upper_gap_RT, ...
    lower_gap_RT, li_RT, dli_dt_RT, q0_RT, qstar_RT, q95_RT, ...
    V_loop_efit_RT, Wmhd_RT, dWmhd_dt_RT] = ...
  get_EFIT_parameters_RT(shot, timebase);

% This function obtains the real-time EFIT parameters from the EFITRT1 tree
% and also calculates the time derivatives of a few of them.
%
% Inputs:
%   shot = shot number
%   timebase = array of times at which to evaluate the parameters
%
% Outputs:
%   beta_N_RT = normalized beta
%   beta_p_RT = beta_poloidal
%   dbeta_p_dt_RT = time derivative of beta_p [/s]
%   kappa_RT = elongation
%   upper_gap_RT = upper gap [m]
%   lower_gap_RT = lower gap [m]
%   li_RT = internal inductance
%   dli_dt_RT = time derivative of li [/s]
%   q0_RT = safety factor at center
%   qstar_RT = q* = cylindrical safety factor
%   q95_RT = edge safety factor
%   V_loop_efit_RT = loop voltage [V]
%   Wmhd_RT = thermal plasma energy, [J]
%   dWmhd_dt_RT = time derivative of Wmhd, [J/s]
%
% Author: Alex Tinguely 2015-10-07 (Dan Brunner's 30th Birthday!)
% Some modifications by Robert Granetz 2015/10/20
% Updates: DIII-D by Alex Tinguely 2016-02-26
% Updates: DIII-D EFIT parameters changes by Cristina Rea 2016-08-11
% Updates: DIII-D real-time EFITRT1 parameters changes by Kevin Montes
%          2017-08-31

% Initialize all outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

beta_N_RT = NaN(size(timebase));
beta_p_RT = NaN(size(timebase));
dbeta_p_dt_RT = NaN(size(timebase));
kappa_RT = NaN(size(timebase));
upper_gap_RT = NaN(size(timebase));
lower_gap_RT = NaN(size(timebase));
li_RT = NaN(size(timebase));
dli_dt_RT = NaN(size(timebase));
q0_RT = NaN(size(timebase));
qstar_RT = NaN(size(timebase));
q95_RT = NaN(size(timebase));
V_loop_efit_RT = NaN(size(timebase));
Wmhd_RT = NaN(size(timebase));
dWmhd_dt_RT = NaN(size(timebase));

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

% For the real-time (RT) signals, read from the EFITRT1 tree; return if the
% tree fails to open.

[~, status] = mdsopen('efitrt1', shot);
if (mod(status,2) == 0);
  return;
end;

% Read in EFIT timebase.  DIII-D is "stupid" (Bob quote) and records time
% in ms. Also, note the efit timebase data is in a node called "atime"
% instead of "time" (where "time" does not work).

[efittime, status] = mdsvalue('\efit_a_eqdsk:atime');
efittime = efittime/1000; % efit time in seconds

if (mod(status,2)==0 || length(efittime) <= 4);
  return;
end;

beta_N_RT = mdsvalue('\efit_a_eqdsk:betan');
beta_p_RT = mdsvalue('\efit_a_eqdsk:betap');
kappa_RT = mdsvalue('\efit_a_eqdsk:kappa'); % changed from "eout" from C-Mod
li_RT = mdsvalue('\efit_a_eqdsk:li');
upper_gap_RT = mdsvalue('\efit_a_eqdsk:gaptop'); % meters, changed from "otop"
                                             % from C-Mod 
lower_gap_RT = mdsvalue('\efit_a_eqdsk:gapbot'); % meters, changed from
                                              % "obott" from C-Mod 
q0_RT = mdsvalue('\efit_a_eqdsk:q0');
qstar_RT = mdsvalue('\efit_a_eqdsk:qstar');
q95_RT = mdsvalue('\efit_a_eqdsk:q95');
V_loop_efit_RT = mdsvalue('\efit_a_eqdsk:vsurf'); % Carlos Paz-Soldan
                                               % recommended that we use the
					       % vsurf signal
Wmhd_RT = mdsvalue('\efit_a_eqdsk:wmhd'); % wmhd is the EFIT reconstruction
                                        % of stored energy 
chisq = mdsvalue('\efit_a_eqdsk:chisq'); % Use chisq to determine which time
                                         % slices are invalid.
mdsclose;

% EFIT reconstructions are sometimes invalid, particularly when very close
% to a disruption.  There are a number of EFIT parameters that can indicate
% invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
% 'chisq' to determine which time slices should be excluded from our
% disruption warning database.

invalid_indices = find(chisq > 50);

beta_N_RT(invalid_indices) = NaN;
beta_p_RT(invalid_indices) = NaN;
kappa_RT(invalid_indices) = NaN;
li_RT(invalid_indices) = NaN;
upper_gap_RT(invalid_indices) = NaN;
lower_gap_RT(invalid_indices) = NaN;
q0_RT(invalid_indices) = NaN;
qstar_RT(invalid_indices) = NaN;
q95_RT(invalid_indices) = NaN;
V_loop_efit_RT(invalid_indices) = NaN;
Wmhd_RT(invalid_indices) = NaN;

dbeta_p_dt_RT = gradient(beta_p_RT, efittime);
dli_dt_RT = gradient(li_RT, efittime);
dWmhd_dt_RT = gradient(Wmhd_RT, efittime);

% Interpolate data onto the requested timebase

beta_N_RT = interp1(efittime, beta_N_RT, timebase_column, 'linear');
beta_p_RT = interp1(efittime, beta_p_RT, timebase_column, 'linear');
kappa_RT = interp1(efittime, kappa_RT, timebase_column, 'linear');
li_RT = interp1(efittime, li_RT, timebase_column, 'linear');
upper_gap_RT = interp1(efittime, upper_gap_RT, timebase_column, 'linear');
lower_gap_RT = interp1(efittime, lower_gap_RT, timebase_column, 'linear');
q0_RT= interp1(efittime, q0_RT, timebase_column, 'linear');
qstar_RT = interp1(efittime, qstar_RT, timebase_column, 'linear');
q95_RT = interp1(efittime, q95_RT, timebase_column, 'linear');
V_loop_efit_RT = interp1(efittime, V_loop_efit_RT, timebase_column, 'linear');
Wmhd_RT = interp1(efittime, Wmhd_RT, timebase_column, 'linear');
dbeta_p_dt_RT = interp1(efittime, dbeta_p_dt_RT, timebase_column, 'linear');
dli_dt_RT = interp1(efittime, dli_dt_RT, timebase_column, 'linear');
dWmhd_dt_RT = interp1(efittime, dWmhd_dt_RT, timebase_column, 'linear');

if (size(timebase,2) > 1);    % If "timebase" is a row vector, then convert
  beta_N_RT = transpose(beta_N_RT); % all the outputs back to row vectors
  beta_p_RT = transpose(beta_p_RT);
  kappa_RT = transpose(kappa_RT);
  li_RT = transpose(li_RT);
  upper_gap_RT = transpose(upper_gap_RT);
  lower_gap_RT = transpose(lower_gap_RT);
  q0_RT = transpose(q0_RT);
  qstar_RT = transpose(qstar_RT);
  q95_RT = transpose(q95_RT);
  V_loop_efit_RT = transpose(V_loop_efit_RT);
  Wmhd_RT = transpose(Wmhd_RT);
  dbeta_p_dt_RT = transpose(dbeta_p_dt_RT);
  dli_dt_RT = transpose(dli_dt_RT);
  dWmhd_dt_RT = transpose(dWmhd_dt_RT);
end;

end
