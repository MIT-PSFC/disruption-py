function [ne, Greenwald_fraction, dne_dt] = ...
  get_density_parameters(shot, timebase);

% This routine obtains the line-averaged density MDS nodename \DENSITY
% (actually a TDI function) in the EFIT trees created for the disruption
% warning database.  It is calculated using the double-pass interferometer
% line-integrated density measurements and EFIT path lengths (PATHV1,
% PATHV2, PATHV3, PATHR0)
% The routine also calculates the Greenwald density using aminor coming
% from our EFIT trees (run_by = {'granetzr', 'reac'}, and runtag = 'DIS').
% Time derivative of electron density is also obtained.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ne = density as measured by the interferometer [m-3]
%   Greenwald_fraction = Greenwald density [m-3]
%   dne_dt = d(ne)/dt [m-3/s]
%
%
% Authors: Cristina Rea and Robert Granetz   Dec 2016
% Revision/history:
%
%  2018/07/09 - RSG; changed to read the data from the EFIT trees that we
%                    created specifically for the disruption warning
%                    database, instead of the default EFIT01 trees (because
%                    we recently discovered that the \density TDI function
%                    stopped working in the EFIT01 runs after 2016/04/22).
%                    Also, if the \density node in our EFIT trees throws an
%                    error, this routine will now try to read the data from
%                    the \denv2 node in the BCI subtree of the D3D main
%                    tree.

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

ne = NaN(size(timebase));
Greenwald_fraction = NaN(size(timebase));
dne_dt = NaN(size(timebase));

mdsconnect('atlas.gat.com');

efittrees = select_efit_trees(shot, '', 'DIS');
if isempty(efittrees);
  efit_disruption_tree_exists = 0;
else;
  tree = char(efittrees(end,:));
  [~, status] = mdsopen(tree, shot);
  if (mod(status,2) == 0);
    efit_disruption_tree_exists = 0;
  else;
    efit_disruption_tree_exists = 1;
  end;
end;

% Read in the density.  If successful, calculate its time derivative.  Then
% interpolate both signals onto the requested timebase.

if efit_disruption_tree_exists;
  [netime, status] = mdsvalue('dim_of(\density)');
  if (mod(status,2) == 1 && length(netime) >= 4);
    netime = netime/1.e3; % convert ms to s
    ne = mdsvalue('\density') * 1.e6; %[m-3]
    mdsclose;
    dne_dt = gradient(ne, netime); %[m-3/s]
    ne = interp1(netime, ne, timebase_column, 'linear');
    dne_dt = interp1(netime, dne_dt, timebase_column, 'linear');
  end;
end;

% If EFIT disruption tree does not exist, or if it does not contain density
% data, then read density from BCI subtree of D3D main tree

if (~efit_disruption_tree_exists || length(find(~isnan(ne))) == 0);
  [~, status] = mdsopen('d3d', shot);
  if (mod(status,2) == 0);
%   fprintf(1,'  Unable to get density data for shot%7i\n', shot);
    return;
  else;
    [netime, status] = mdsvalue('dim_of(\denv2)');
    if (mod(status,2) == 0);
%     fprintf(1,'  Unable to get density data for shot%7i\n', shot);
      return;
    else;
      netime = netime/1.e3; % convert ms to s
      ne = mdsvalue('\denv2') / 2.774 * 1.e6; % convert to m-3
      mdsclose;
      dne_dt = gradient(ne, netime); %[m-3/s]
      ne = interp1(netime, ne, timebase_column, 'linear');
      dne_dt = interp1(netime, dne_dt, timebase_column, 'linear');
    end;
  end;
end;

% Next, get the plasma current, which is needed to calculate the Greenwald
% density

[iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' num2str(shot) '))']);% ms
if (mod(status,2) == 1);
  iptime = iptime/1.e3; % convert ms to s
  ip = mdsvalue(['ptdata("ip", ' num2str(shot) ')']); %[A]
  ipsign = sign(sum(ip));
  ip = interp1(iptime, ip*ipsign, timebase_column, 'linear');%positive definite
else;
  if (size(timebase,2) > 1);
    ne = transpose(ne);
    dne_dt = transpose(dne_dt);
  end;
  return;
end;

mdsclose;

efittrees = select_efit_trees(shot,'granetzr','DIS');
if isempty(efittrees);
% fprintf(1, 'No disruption EFIT tree for this shot\n');
  [~, status] = mdsopen('efit01', shot);
  if (mod(status,2) == 0);
    return;
  end;
else;
  tree=char(efittrees(end,:));
  [~, status] = mdsopen(tree, shot);
  if (mod(status,2) == 0);
    return;
  end;
end;

% Read in EFIT minor radius and timebase.  DIII-D records time in ms.
% Also, note the efit timebase data is in a node called "atime" instead of
% "time" (where "time" does not work).

[efittime, status] = mdsvalue('\efit_a_eqdsk:atime');
efittime = efittime/1.e3; % efit time in seconds

if (mod(status,2)==0 || length(efittime) <= 4);
  mdsclose;
  return;
end;

aminor = mdsvalue('\efit_a_eqdsk:aminor'); % meters

mdsclose;

% Interpolate data onto the requested timebase

aminor = interp1(efittime, aminor, timebase_column, 'linear');

nG = (ip/1.e6) ./ (pi*aminor.^2);  %[MA/m2] -> Greenwald ne limit [10^20 m-3]
Greenwald_fraction = (ne/1.e20) ./ nG;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  ne = transpose(ne);
  Greenwald_fraction = transpose(Greenwald_fraction);
  dne_dt = transpose(dne_dt);
end;

end
