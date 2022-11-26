function [ne, nG, dn_dt] = ...
  get_density_parameters(shot, timebase);

% This routine obtains the line-averaged density pointname DENSITY 
% (actually a TDI function), calculated using the double-pass 
% interferometer line-integrated density measurements and EFIT01 path 
% lengths (PATHV1, PATHV2, PATHV3, PATHR0). 
% The routine also calculates the Greenwald density using aminor coming
% from EFIT trees run by Bob ('dis' tagged). 
% Time derivative of electron density is also obtained.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ne = density as measured by the interferometer [m-3]
%   nG = Greenwald density [10^20 m-3]
%   dn_dt = d(ne)/dt [m-3/s]
%
%
% Authors: Cristina Rea and Robert Granetz   Dec 2016

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

ne         = NaN(length(timebase), 1);
nG         = NaN(length(timebase), 1);
dn_dt      = NaN(length(timebase), 1);

mdsconnect('atlas.gat.com');

[shotopened, status]=mdsopen('d3d', shot);
if (mod(status,2)==0);
% fprintf(1,'  Unable to open D3D tree for shot%7i\n', shot);
  return;
end;

% Next, read in the measured plasma current, Ip.
% If Ip is successfully read in, then calculate dI/dt, and then interpolate
% both signals onto the requested timebase.

ne = mdsvalue(['ptdata("density", ' num2str(shot) ')']);%[cm-3]
ne = ne*1.e6; %[m-3]
[netime, status] = mdsvalue(['dim_of(ptdata("density", ' num2str(shot) '))']);
if (mod(status,2) == 1);
  netime = netime/1.e3; % convert ms to s
  dne_dt = gradient(ne, netime); %[m-3/s]
  ne = interp1(netime, ne, timebase_column, 'linear');
  dne_dt = interp1(netime, dne_dt, timebase_column, 'linear');
else;
  ne = NaN(length(timebase), 1);
  dne_dt = NaN(length(timebase), 1);
end;

% Next, get the plasma current to obtain the Greenwald density

ip = mdsvalue(['ptdata("ip", ' num2str(shot) ')']); %[A]
[iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' ...
                                                   num2str(shot) '))']);
if (mod(status,2) == 1);
  iptime = iptime/1.e3; % convert ms to s
  polarity = unique(mdsvalue(['ptdata("iptdirect", ' ...
                                                   num2str(shot) ')']));
  if (length(polarity) ~= 1);
    fprintf(1,'ERROR: polarity of Ip target is not constant\n');
    polarity = polarity(1);
  end;
  ip = ip * polarity;
  ip = interp1(iptime, ip, timebase_column, 'linear');
else;
  ip = NaN(length(timebase), 1);
end;

mdsclose;

efittrees = select_efit_trees(shot,'granetzr','DIS');
if isempty(efittrees);
  fprintf(1, 'No disruption EFIT tree for this shot\n');
  return;
end;
tree=char(efittrees(end,:));
[shotopened, status] = mdsopen(tree, shot); % changed from "efit18" for DIII-D
if (mod(status,2) == 0);
  return;
end;

% Read in EFIT timebase.  DIII-D is "stupid" (Bob quote) and records time
% in ms. Also, note the efit timebase data is in a node called "atime"
% instead of "time" (where "time" does not work).

[efittime, status] = mdsvalue('\efit_a_eqdsk:atime');
efittime = efittime/1.e3; % efit time in seconds

if (mod(status,2)==0 || length(efittime) <= 4);
  return;
end;

aminor = mdsvalue('\efit_a_eqdsk:aminor');

mdsclose;

% Interpolate data onto the requested timebase

aminor = interp1(efittime, aminor, timebase_column, 'linear');

nG = (ip/1.e6) / (pi*(aminor)^2) %[MA/m] -> nG [10^20 m-3]

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  ne         = transpose(ne);
  nG         = transpose(nG);
  dne_dt     = transpose(dne_dt);
end;

end
