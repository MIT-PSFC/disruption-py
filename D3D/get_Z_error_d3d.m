function [z_error, z_error_norm, z_prog, zcur, zcur_norm] = ...
    get_Z_error_d3d(shot, timebase);

% This script was adapted from get_Z_parameters.m on C-Mod.
%
% The purpose of this script is to read in the values of Z_error and Z_prog
% from the plasma control system (PCS). Z_prog is the programmed vertical
% position of the plasma current centroid, and Z_error is the difference
% between the actual position and that requested (Z_error = Z_cur -
% Z_prog). Thus, the actual (estimated) position, Z_cur, can be calculated.
% We do not yet know the programmed Z_prog.
% These values are then linearly interpolated over the given timebase.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values [s]
%
% Outputs:
%   z_error = Z_cur - Z_programmed (m); returns NaN array if no Z_cur or Z_prog
%   z_prog = programmed/requested/target Z (m); returns NaN array if no data
%   zcur = calculated Z from PCS [m]; returns NaN array if no data
%   zcur_norm = zcur / a
%   z_error_norm = z_error / a
%
% Author: Updated by Cristina Rea 2019/08/23 for DIII-D.
%         The routine will return only zcur and zcur_norm properly populated.
%         We still don't know z_prog pointname.
%
% Revision history
%   2019/09/04 - RSG; We found out from Jayson Barr that we misunderstood
%                the PCS pointnames.  The correct pointname is "vpszp".
%   2019/09/23 - RSG; On DIII-D the plasma control system uses isoflux
%                control to control the plasma shape and position.  It does
%                NOT use zcur control.  Therefore, the PCS does not have a
%                programmed vertical position.  This this routine will now
%                always return an arrays of NaN for z_prog, z_error, and
%                z_error_norm.

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

zcur = NaN(length(timebase),1);
zcur_norm = NaN(length(timebase),1);
z_prog = NaN(length(timebase),1);
z_error = NaN(length(timebase),1);
z_error_norm = NaN(length(timebase),1);
aminor = NaN(length(timebase_column),1);

mdsconnect('atlas.gat.com');

[shotopened, status]=mdsopen('d3d', shot);
if (mod(status,2)==0);
% fprintf(1,'  Unable to open D3D tree for shot%7i\n', shot);
  if (size(timebase,2) > 1);
    zcur         = transpose(zcur);
    z_prog    = transpose(z_prog);
    z_error   = transpose(z_error);
    zcur_norm     = transpose(zcur_norm);
    z_error_norm = transpose(z_error_norm);
  end;
  return;
end;

% Next, read in the measured vertical position of the current centroid (PCS).
zcur = mdsvalue(['ptdata("vpszp", ' num2str(shot) ')']);
[ztime, status] = mdsvalue(['dim_of(ptdata("vpszp", ' num2str(shot) '))']);
if (mod(status,2) == 1);
  ztime = ztime/1.e3; % convert ms to s
  zcur = zcur/1.e2; % convert cm to m

% Interpolate data onto the requested timebase
  zcur = interp1(ztime, zcur, timebase_column, 'linear');
else;
  zcur = NaN(length(timebase), 1);
end;

mdsclose;

efittrees  = select_efit_trees(shot, '', 'DIS');
if isempty(efittrees);
  fprintf(1, 'No disruption EFIT tree for this shot\n');
  return;
end;
tree=char(efittrees(end,:));
[shotopened, status] = mdsopen(tree, shot);
if (mod(status,2) == 0);
  return;
end;

% Read in EFIT timebase.  DIII-D is "stupid" (Bob quote) and records time
% in ms. Also, note the efit timebase data is in a node called "atime"
% instead of "time" (where "time" does not work).

[efittime, status] = mdsvalue('\efit_a_eqdsk:atime');

if (mod(status,2)==1 && length(efittime) > 4);
  efittime = efittime/1.e3; % convert ms to s
  aminor = mdsvalue('\efit_a_eqdsk:aminor'); % meters
  chisq = mdsvalue('\efit_a_eqdsk:chisq'); % Use chisq to determine which time
                                           % slices are invalid.

% EFIT reconstructions are sometimes invalid, particularly when very close
% to a disruption.  There are a number of EFIT parameters that can indicate
% invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
% 'chisq' to determine which time slices should be excluded from our
% disruption warning database.

  invalid_indices = find(chisq > 50);
  aminor(invalid_indices) = NaN;
  mdsclose;

% Interpolate data onto the requested timebase and calculate zcur_norm
  aminor = interp1(efittime, aminor, timebase_column, 'linear');
  zcur_norm = zcur ./ aminor;

else;
  mdsclose;
  zcur_norm = zcur / 0.59;
end;

if (size(timebase,2) > 1);    % If "timebase" is a row vector, then convert
  zcur = transpose(zcur); % all the outputs back to row vectors
  zcur_norm = transpose(zcur_norm);
  z_prog = transpose(z_prog);
  z_error = transpose(z_error);
  z_error_norm = transpose(z_error_norm);
end;

end
