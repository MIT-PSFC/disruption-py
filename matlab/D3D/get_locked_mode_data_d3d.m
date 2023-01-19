function [n_equal_1_mode, IRLM] = get_locked_mode_data_d3d(shot, timebase);

% The original locked mode data from Ryan Sweeney is written in an IDL
% saveset file named "irlms_dataarray.sav".  The publicly-available Matlab
% program "restore_idl.m" can read a select set of IDL variable types from
% IDL saveset files and convert them to Matlab structures.  HOWEVER, Ryan's
% saveset contains pointer heap variables, which cannot be handled by the
% "restore_idl.m" Matlab routine.  So first, an IDL program that I wrote,
% named "read_locked_mode_idl_saveset.pro", must be run to convert Ryan's
% IDL saveset into a new saveset, "locked_mode_idl_saveset.sav", which
% contains all the data in variable types that "restore_idl.m" supports.
% After that is done, then this routine, which uses "restore_idl.m" can be
% called.  Note that "restore_idl.m" returns a structure whose field names
% are uppercase by default.

persistent SHOTS PTS TIMES LOCKED_MODE_AMPLITUDE

if isempty(SHOTS);
  outstruct = restore_idl(['/fusion/projects/disruption_warning/' ...
    'idl_programs/locked_mode_idl_saveset.sav']);

% structure fields: SHOTS, PTS, TIMES, LOCKED_MODE_AMPLITUDE

  SHOTS = outstruct.SHOTS;
  PTS = outstruct.PTS;
  TIMES = outstruct.TIMES;
  LOCKED_MODE_AMPLITUDE = outstruct.LOCKED_MODE_AMPLITUDE;
end;

[~, indx_shot] = ismember(shot, SHOTS);
if indx_shot ~= 0;
  npts = PTS(indx_shot);
  time = TIMES(indx_shot, 1:npts);
  n_equal_1_mode = LOCKED_MODE_AMPLITUDE(indx_shot, 1:npts);

% So now we have the n=1 locked mode amplitude vs time for this shot.  But
% these data may not be continuous in time.  They may be separated into
% blocks of time, with no information available during the times between
% blocks.  However, we want to interpolate values at all of the requested
% times ('times_for_db'), which likely include the gap(s), if any, between
% blocks of data.  Therefore, we are going to fill in any gaps with the
% default value, which is zero.

  dt = min(diff(time));
  jump_indices = find(diff(time) > 1.1*dt);
  if length(jump_indices) > 0;
    for jump = 1:length(jump_indices);
      jump_indx = jump_indices(jump) + 2*(jump-1);
      tjump1 = time(jump_indx);
      tjump2 = time(jump_indx+1);
      time = [time(1 : jump_indx), tjump1+1.e-6, tjump2-1.e-6, ...
              time(jump_indx+1 : end)];
      n_equal_1_mode = [n_equal_1_mode(1 : jump_indx), 0, 0, ...
                        n_equal_1_mode(jump_indx+1 : end)];
    end;
  end;

% Make output arrays the same shape as "timebase" (i.e. row or column)
  if (size(timebase,1) > 1) && (size(time,1) == 1);
    time = transpose(time);
  end;
  if (size(timebase,1) > 1) && (size(n_equal_1_mode,1) == 1);
    n_equal_1_mode = transpose(n_equal_1_mode);
  end;

% Interpolate onto requested timebase
  n_equal_1_mode = interp1(time, n_equal_1_mode, timebase, 'linear', 0.0);

  IRLM = ones(size(timebase));
else;
  n_equal_1_mode = zeros(size(timebase));
  IRLM = zeros(size(timebase));
end;

return;
end
