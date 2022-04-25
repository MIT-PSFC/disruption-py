function [ne, Greenwald_fraction, dn_dt] = ...
  get_density_parameters(shot, timebase);

% This routine obtains the line-averaged density from the HCN vertical
% chord.  The node is called \DFSDEV in the PCS_EAST tree.  This signal is
% also used by the PCS for feedback control of the density.
%
% The plasma density is also measured by a POlarimeter/INTerferometer
% diagnostic (PO INT), which has 11 horizontal chords.  The midplane
% channel is \POINT_N6 in the EAST tree.  (This diagnostic uses an FIR
% laser at 432 micrometers.)  However, I found that the POINT density
% measurement is bad on a significant fraction of EAST shots, even in
% 2017.
%
% The routine also calculates the Greenwald density using aminor coming
% from EFIT.  The time derivative of electron density is also obtained.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ne = density as measured by the polaritmeter/interferometer [m-3]
%   Greenwald_fraction = ne/nG
%   dn_dt = d(ne)/dt [m-3/s]
%
% Author: Robert Granetz   Dec 2017

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

ne                 = NaN(length(timebase), 1);
Greenwald_fraction = NaN(length(timebase), 1);
dn_dt              = NaN(length(timebase), 1);

mdsconnect('mds.ipp.ac.cn');

% Get the density

[shotopened, status] = mdsopen('pcs_east', shot);
if (mod(status,2) == 0);
% fprintf(1,'  Unable to open EAST tree for shot%7i\n', shot);
  return;
end;

% Read in the measured density.  If ne is successfully read in, then
% calculate dn/dt, and then interpolate both signals onto the requested
% timebase.

[netime, status] = mdsvalue('dim_of(\dfsdev)');
if (mod(status,2) == 1);
  ne = mdsvalue('\dfsdev'); % density [10^19 m-3]
  ne = ne * 1.e19; % [m-3]
else;
  return;
end;

% Don't use the POINT_N6 data
%{

[shotopened, status]=mdsopen('east', shot);
if (mod(status,2) == 0);
% fprintf(1,'  Unable to open EAST tree for shot%7i\n', shot);
  return;
end;

[netime, status] = mdsvalue('dim_of(\point_n6)');
if (mod(status,2) == 1);
  neL = mdsvalue('\point_n6'); % nL [10^19 m-2]
  neL = neL*1.e19; % [m-2]
else;
  return;
end;
%}

% I need the minor radius, aminor, to calculate the Greenwald density.  I
% get this from EFIT (\efit_aeqdsk:aout)

% Open EFIT18 tree.  If error, try opening the standard EFIT tree,
% 'EFIT_EAST'.  If that also fails, then use a default value of 0.45 m

[shotopened, status] = mdsopen('efit18', double(shot));
if (mod(status,2) == 1);
  [efittime, status] = mdsvalue('dim_of(\efit_aeqdsk:aout)');
  if (mod(status,2) == 1 & length(efittime) >= 2);
    [aminor, status] = mdsvalue('\efit_aeqdsk:aout');
    if (mod(status,2) == 0);
      aminor = 0.45 * ones(size(efittime));
    end;
  else;
    efittime = [0; 1000]; % If no efit timebase, make fake one
    aminor = [0.45; 0.45];
  end;
  mdsclose;
else;
  [shotopened, status] = mdsopen('efit_east', double(shot));
  if (mod(status,2) == 1);
    [efittime, status] = mdsvalue('dim_of(\efit_aeqdsk:aout)');
    if (mod(status,2) == 1 & length(efittime) >= 2);
      aminor = mdsvalue('\efit_aeqdsk:aout');
    else;
      efittime = [0; 1000]; % If no efit timebase, make fake one
      aminor = [0.45; 0.45];
    end;
    mdsclose;
  else;
    efittime = [0; 1000]; % If no efit timebase, make fake one
    aminor = [0.45; 0.45];
  end;
end;

%{
if (length(efittime) > length(aminor));    % bug: efittime sometimes has
  efittime = efittime(1 : length(aminor)); % an extra element (KSTAR)
end;
%}

[efittime, indices_unique] = unique(efittime); % Another rare bug: sometimes
aminor = aminor(indices_unique);               % the efit timebase has a
                                               % repeated value

%{
% Don't use the POINT_N6 data

% Need to divide nL by L to get <n>.  But first, I must interpolate aminor
% onto the neL timebase

aminor_on_neL_timebase = interp1(efittime, aminor, netime, 'linear');

% Now divide nL by L to get <n>

ne = neL ./ (2*aminor_on_neL_timebase);
%}

dn_dt = gradient(ne, netime); % [m-3/s]

ne = interp1(netime, ne, timebase_column, 'linear');
dn_dt = interp1(netime, dn_dt, timebase_column, 'linear');

% Next, get the plasma current so I can calculate the Greenwald density.
% I also need the plasma minor radius, but I already have that.

% Get \PCRL01 plasma current signal from the PCS_EAST tree.  If there is no
% PCS_EAST tree, or no PCRL01 data, then this is not a valid plasma shot
% (as per guidance from Qian Jinping).

[shotopened, status] = mdsopen('pcs_east', double(shot));
if (mod(status,2) == 1);
  [ip, status] = mdsvalue('\pcrl01');  % Read in Ip data (amps)
  if (mod(status,2) == 1);             % If successful, continue
    iptime = mdsvalue('dim_of(\pcrl01)');

% For shots before year 2014, the PCRL01 timebase needs to be shifted by
% 17.0 ms 

    if (shot < 44432);
      iptime = iptime - 0.0170;
    end;

    baseindices=find(iptime <= -5.8); % time before any PF supplies turn on
    if (numel(baseindices) > 0);
      baseline=sum(ip(baseindices))/length(baseindices);
      ip = ip - baseline;
    end;

%{
    if (length(iptime) > length(ip));  % bug: timebase occasionally
      iptime = iptime(1 : length(ip)); % has an extra element (KSTAR)
    end;
%}

    ip = interp1(iptime, ip, timebase_column, 'linear');
    mdsclose;
  else;
    if (size(timebase,2) > 1);
      ne    = transpose(ne);
      dn_dt = transpose(dn_dt);
    end;
    return;
  end;
else;
  if (size(timebase,2) > 1);
    ne    = transpose(ne);
    dn_dt = transpose(dn_dt);
  end;
  return;
end;

aminor = interp1(efittime, aminor, timebase_column, 'linear');

% Now calculate Greenwald density, then Greenwald fraction.  But first,
% if the plasma current is negative (normal direction), then invert it,
% since the empirical formula for the Greenwald density doesn't depend on
% the sign of the plasma current.

if (sum(ip) < 0); ip = -ip; end;

nG = (ip/1.e6) ./ (pi*aminor.^2); % [MA/m2] -> nG in 10^20 m-3
nG = nG * 1.e20; % convert to m-3

Greenwald_fraction = ne ./ nG;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  ne                 = transpose(ne);
  Greenwald_fraction = transpose(Greenwald_fraction);
  dn_dt              = transpose(dn_dt);
end;

end
