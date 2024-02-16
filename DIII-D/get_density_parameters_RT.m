function [ne_RT, Greenwald_fraction_RT, dne_dt_RT] = ...
  get_density_parameters_RT(shot, timebase);

% This routine obtains the line-averaged density available to the PCS in
% real time ('dssdenest'), which is calculated using the double-pass
% interferometer line-integrated density measurements and EFITRT1 path
% lengths (PATHV1, PATHV2, PATHV3, PATHR0).  This routine also calculates
% the Greenwald density using aminor coming from EFITRT1 trees.  The time
% derivative of electron density is also obtained.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ne_RT = density as measured by the interferometer [m-3]
%   Greenwald_fraction_RT = Greenwald density [m-3]
%   dne_dt_RT = d(ne)/dt [m-3/s]
%
%
% Authors: Cristina Rea and Robert Granetz   Aug 2017

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

ne_RT = NaN(size(timebase));
Greenwald_fraction_RT = NaN(size(timebase));
dne_dt_RT = NaN(size(timebase));

% Read in the PCS real-time density signal.  If successful, calculate its
% time derivative.  Then interpolate both signals onto the requested
% timebase.

[netime, status] = mdsvalue(['dim_of(ptdata("dssdenest", ' ...
  num2str(shot) '))']);
if (mod(status,2) == 1);
  netime = netime/1.e3; % convert ms to s
  ne_RT =  mdsvalue(['ptdata("dssdenest", ' num2str(shot) ')']); %[10^19 m-3]
  ne_RT = ne_RT*1.e19; %[m-3]
  dne_dt_RT = gradient(ne_RT, netime); %[m-3/s]
  ne_RT = interp1(netime, ne_RT, timebase_column, 'linear');
  dne_dt_RT = interp1(netime, dne_dt_RT, timebase_column, 'linear');
else;
  return;
end;

% Next, get the plasma current, which is needed to calculate the Greenwald
% density

[iptime, status]=mdsvalue(['dim_of(ptdata("ipsip", ' num2str(shot) '))']);% ms
if (mod(status,2) == 1);
  iptime = iptime/1.e3; % convert ms to s
  [ip, status] = mdsvalue(['ptdata("ipsip", ' num2str(shot) ')']); %[MA]
else;
  [iptime, status] = mdsvalue(['dim_of(ptdata("ipspr15v" ', ...
    num2str(shot) '))']);% ms
  if (mod(status,2) == 1);
    iptime = iptime/1.e3; % convert ms to s
    ip = mdsvalue(['ptdata("ipspr15v", ' num2str(shot) ')']); %[volts; 2 V/MA]
    ip = ip/2; % Convert volts to MA
  else;
    return;
  end;
end;

ipsign = sign(sum(ip));
ip = interp1(iptime, ip*ipsign, timebase_column, 'linear');%positive definite

% Read in EFIT minor radius and timebase.  This is also needed to calculate
% the Greenwald density limit.  However, if the minor radius data is not
% available, use a default fixed value of 0.59 m.  (We surveyed several
% hundred shots to determine this default value.)  Note that the efit
% timebase data is in a node called "atime" instead of "time" (where "time"
% does not work).

% For the real-time (RT) signals, read from the EFITRT1 tree

[~, status] = mdsopen('efitrt1', shot);
if (mod(status,2) == 1);
  [efittime, status] = mdsvalue('\efit_a_eqdsk:atime');
  if (mod(status,2) == 1);
    efittime = efittime/1.e3; % efit time in seconds
    aminor = mdsvalue('\efit_a_eqdsk:aminor'); % meters
% Interpolate data onto the requested timebase
    aminor = interp1(efittime, aminor, timebase_column, 'linear');
  else;
    aminor = 0.59 * ones(size(timebase_column));
  end;
  mdsclose;
else;
  aminor = 0.59 * ones(size(timebase_column));
end;

nG = ip ./ (pi*aminor.^2);  %[MA/m2] -> Greenwald ne limit [10^20 m-3]
Greenwald_fraction_RT = (ne_RT/1.e20) ./ nG;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  ne_RT = transpose(ne_RT);
  Greenwald_fraction_RT = transpose(Greenwald_fraction_RT);
  dne_dt_RT = transpose(dne_dt_RT);
end;

end
