function [Z_error, Z_prog, Z_cur, zcur_lmsz, z_error_lmsz, ...
  zcur_lmsz_normalized, z_error_lmsz_normalized] = ...
  get_Z_error_EAST(shot, timebase);

% This script calculates Z_error = Z_cur - Z_programmed, or how much the
% actual vertical position differs from the requested position.  Two
% different methods are used to calculate Z_error, and both versions are
% returned by the routine.  The original method gets the z-centroid from
% the \zcur calculated by EFIT; the other uses the z-centroid from the
% \lmsz signal calculated by the plasma control system (PCS).  The
% normalized versions of the lmsz-derived signals are also returned.  They
% are normalized to the plasma minor radius.  All the signals are linearly
% interpolated onto the given timebase.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   Z_error = Z_cur - Z_programmed (m); returns NaN array if no Z_cur or Z_prog
%   Z_prog = programmed/requested/target Z (m); returns NaN array if no data
%   Z_cur = calculated Z from EFIT [m]; returns NaN array if no data
%   zcur_lmsz = calculated Z from PCS [m]; returns NaN array if no data
%   z_error_lmsz = zcur_lmsz - Z_programmed [m]; returns NaN array if no data
%   zcur_lmsz_normalized = zcur_lmsz / a
%   z_error_lmsz_normalized = z_error_lmsz / a
%
% Author: Modified on EAST by Wang Bo 2015/12/10
%         Alex Tinguely 2015-09-09
%         Some modifications by RS Granetz 
%
% Revision history:
%  2018/08/08 - RSG; Added output signals based on \lmsz from PCS.
%                    Also modified program to read \zcur from our EFIT18
%                    tree, if it exists.
%  2019/05/10 - RSG; Added output of normalized versions of the lmsz-based
%                    signals.  They are normalized by dividing by the
%                    plasma minor radius, which is obtained from EFIT.

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

% Initialize all output arrays to NaN (Not-a-Number) column vectors

Z_cur   = NaN(size(timebase_column));
Z_prog  = NaN(size(timebase_column));
Z_error = NaN(size(timebase_column));
zcur_lmsz = NaN(size(timebase_column));
z_error_lmsz = NaN(size(timebase_column));
zcur_lmsz_normalized = NaN(size(timebase_column));
z_error_lmsz_normalized = NaN(size(timebase_column));
aminor = NaN(size(timebase_column));

% Next, read in the calculated Zcur from EFIT.  Try our EFIT18 tree first,
% and if that fails, use the standard EFIT_EAST tree.

[~, status] = mdsopen('efit18', double(shot));
if (mod(status,2) == 1);
  [Z_cur_time, status] = mdsvalue('dim_of(\efit_aeqdsk:zcur)');
  if (mod(status,2) == 1 && length(Z_cur_time) > 2);     % deal with rare bug
    [Z_cur_time, unique_indices, ~] = unique(Z_cur_time);% deal with rare bug
    Z_cur = mdsvalue('\efit_aeqdsk:zcur');  % meters
    Z_cur = Z_cur(unique_indices);                       % deal with rare bug
% Interpolate onto requested timebase
    Z_cur = interp1(Z_cur_time, Z_cur, timebase_column, 'linear');

    [aminor, status] = mdsvalue('\aminor');  % needed for normalization
    if (mod(status,2) == 1);
      aminor = aminor(unique_indices);
      aminor = interp1(Z_cur_time, aminor, timebase_column, 'linear');
    end;
  end;
  mdsclose;
else;   % If unable to read from our EFIT18 tree, then:
  [~, status] = mdsopen('efit_east', double(shot));
  if (mod(status,2) == 1);
    [Z_cur_time, status] = mdsvalue('dim_of(\efit_aeqdsk:zcur)');
    if (mod(status,2) == 1 && length(Z_cur_time) > 2);   % deal with rare bug
      [Z_cur_time, unique_indices, ~] = unique(Z_cur_time); % deal w/rare bug
      Z_cur = mdsvalue('\efit_aeqdsk:zcur');  % meters
      Z_cur = Z_cur(unique_indices);                     % deal with rare bug
% Interpolate onto requested timebase
      Z_cur = interp1(Z_cur_time, Z_cur, timebase_column, 'linear');

      [aminor, status] = mdsvalue('\aminor');  % needed for normalization
      if (mod(status,2) == 1);
        aminor = aminor(unique_indices);
        aminor = interp1(Z_cur_time, aminor, timebase_column, 'linear');
      end;
    end;
    mdsclose;
  end;
end;

% Next, get the programmed/requested/target Z from PCS, and the
% calculated Z-centroid from PCS

[shotopened,status] = mdsopen('pcs_east', double(shot));
if (mod(status,2) == 1);
  [Z_prog_time, status] = mdsvalue('dim_of(\lmtzref)');
  if (mod(status,2)==1 && length(Z_prog_time) > 2);
    Z_prog = mdsvalue('\lmtzref');  % cm (node says 'm', but it's wrong)
    Z_prog = Z_prog/100; % convert to meters
% Interpolate onto requested timebase
    Z_prog = interp1(Z_prog_time, Z_prog, timebase_column, 'linear');
  end;

  [lmsz_time, status] = mdsvalue('dim_of(\lmsz)');
  if (mod(status,2)==1 && length(lmsz_time) > 2);
    zcur_lmsz = mdsvalue('\lmsz'); % [m]
% Interpolate onto requested timebase
    zcur_lmsz = interp1(lmsz_time, zcur_lmsz, timebase_column, 'linear');
    zcur_lmsz_normalized = zcur_lmsz ./ aminor;
  end;
  mdsclose;
end;

% Calculate both versions of Z_error

Z_error = Z_cur - Z_prog;  % meters
z_error_lmsz = zcur_lmsz - Z_prog;  % meters
z_error_lmsz_normalized = z_error_lmsz ./ aminor;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Z_cur   = transpose(Z_cur);
  Z_prog  = transpose(Z_prog);
  Z_error = transpose(Z_error);
  zcur_lmsz = transpose(zcur_lmsz);
  z_error_lmsz = transpose(z_error_lmsz);
  zcur_lmsz_normalized = transpose(zcur_lmsz_normalized);
  z_error_lmsz_normalized = transpose(z_error_lmsz_normalized);
end;

end
