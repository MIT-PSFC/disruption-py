function sxr = get_sxr_data(shot, timebase);

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

[~, status] = mdsopen('xtomo', shot);
if (mod(status, 2) == 0);
  sxr = NaN(size(timebase));
  return;
end;

[sxr_time, status] = mdsvalue('dim_of(\top.brightnesses.array_1:chord_16)');
if (mod(status, 2) == 1 && length(sxr_time) > 1);
  sxr = mdsvalue('\top.brightnesses.array_1:chord_16');
else;
  sxr = NaN(size(timebase));
  return;
end;
mdsclose;

sxr = interp1(sxr_time, sxr, timebase_column, 'linear'); 

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase, 2) > 1);
  sxr = transpose(sxr);
end;

end