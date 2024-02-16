function [Te_HWHM_RT, Te_width_normalized_RT] = get_TS_data_RT(shot, timebase);

% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Output:
%   Te_HWHM_RT(timebase) = half-width at half-max of parabolic fits to the
%     Thomson Te(t,z) data available in real time from the DIII-D PCS
%   Te_width_normalized_RT(timebase) = HWHM normalized to minor radius
%
% Authors: Robert Granetz and Cristina Rea        Aug 2017
%     Aug 2017 - modified to use the Thomson data as recorded by PCS,
%                instead of the "blessed" TS data stored in the MDS tree
%     Jan 2018 - added normalized width

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

% Read in Thomson core temperature data, which is a 2-D array, with the
% dependent dimensions being time and z (vertical coordinate)

[TS_time, status] = mdsvalue(['dim_of(ptdata("tsscorte1", ' ...
  num2str(shot) '))']);

if (mod(status,2) == 1 & length(TS_time) >= 2);
  TS_data = NaN(length(TS_time), 40);
  for ichan = 1:40;
    TS_single_channel = mdsvalue(['ptdata("tsscorte' num2str(ichan) '", ' ...
      num2str(shot) ')']);
    TS_data(:, ichan) = TS_single_channel;
  end;

% Get the TS z-coordinates.  These are not available to the PCS, but they
% are stored in the MDS tree for each shot.
% Oct 10 2019: Note that R,Z location of Thomson lasers is subject to variations
% year-by-year as well as after mainteneance downtime.
% DO NOT RELY on old shots to provide R, Z coordinates.

  mdsconnect('atlas.gat.com');

  [shotopened, status]=mdsopen('electrons', shot);
  if (mod(status,2)==1);
    [TS_z, status] = mdsvalue('dim_of(\top.ts.blessed.core:temp,1)');
    if (mod(status,2) ~= 1);
      load('TS_z_coordinate_data.mat'); % get saved TS z coords from 164716
    end;
  else;
    fprintf(1,'  Unable to open ELECTRONS tree for shot%7i\n', shot);
    load('TS_z_coordinate_data.mat'); % get saved TS z coordinates from 164716
  end;
  mdsclose;

  TS_time = TS_time/1.e3; % convert ms to s

  Te_HWHM_RT = NaN(length(TS_time),1);

  zarray = [0:.01:.9];
  itimes = find(TS_time > 0);

  for i = 1:length(itimes);
    y = TS_data(itimes(i),:);
    ok_indices = find(y ~= 0);
    y = y(ok_indices);
    z = TS_z(ok_indices);
    if (length(ok_indices) > 2);
      p = polyfit(z, transpose(y), 2);
      Te_array = polyval(p, zarray);
      [Te_max, maxindx] = max(Te_array);
      z_max = zarray(maxindx);
      Te_HM = Te_max/2;
      [~, HM_indices] = min(abs(Te_array - Te_HM));
      HM_indx = max(HM_indices);
      z_HM = zarray(HM_indx);
      if z_HM > z_max;
        Te_HWHM_RT(itimes(i)) = z_HM - z_max;
      end;
    end;
  end;
  Te_HWHM_RT = interp1(TS_time, Te_HWHM_RT, timebase_column, 'linear');
else;
  mdsclose;
  Te_HWHM_RT = NaN(size(timebase));
  Te_width_normalized_RT = NaN(size(timebase));
  return;
end;


% Read in EFIT minor radius and timebase.  This is needed to normalize the
% Te width.  However, if the minor radius data is not available, use a
% default fixed value of 0.59 m.  (We surveyed several hundred shots to
% determine this default value.)  Note that the efit timebase data is in a
% node called "atime" instead of "time" (where "time" does not work).

% For the real-time (RT) signals, read from the EFITRT1 tree

[~, status] = mdsopen('efitrt1', shot);
if (mod(status,2) == 1);
  [efittime, status] = mdsvalue('\efit_a_eqdsk:atime');
  if (mod(status,2) == 1 && length(efittime) > 4);
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

Te_width_normalized_RT = Te_HWHM_RT ./ aminor;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Te_HWHM_RT = transpose(Te_HWHM_RT);
  Te_width_normalized_RT = transpose(Te_width_normalized_RT);
end;

end
