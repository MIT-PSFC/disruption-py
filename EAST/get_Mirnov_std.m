function [Mirnov_std_normalized, Mirnov_std] = ...
  get_Mirnov_std(shot, timebase, Mirnov_sensor)

% the set of Mirnov sensors digitized at 200 kHz (from Dalong):
%
% Name: cmp1t~cmp26t  (c-port), kmp1t~kmp26t (k-port),
% all in the EAST tree

if (nargin<3);
  Mirnov_sensor = 'cmp1t';
end;
 
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

ntimes = length(timebase);
Mirnov_std = NaN(size(timebase));
Mirnov_std_normalized = NaN(size(timebase));

% Next, get the toroidal magnetic field.  For shots < 60000, the TF current
% was in node "\it" in the "eng_tree", and the timebase of the signal was
% not well-defined.  For shots between 60000 (on 2016/01/29) and 65165
% (2016/05/01), the "\it" node is in the "pcs_east" tree, with a proper
% timebase.  For shots after 65165, the "\it" node was back in the
% "eng_tree" tree (with a proper timebase).  Also, starting with shot
% 60000, there is a fibreoptic-based measurement of the TF current.  The
% signals are "\focs_it" (digitized at 50 kHz) and "\focs_it_s"
% (sub-sampled at 1 kHz), both in the "east" tree.  The latter two signals
% differ by 1.6% from the \it signal (as of 2016/04/18).

if (shot > 65165);
  [~, status] = mdsopen('eng_tree', double(shot));
  if (mod(status, 2) == 1);
    [btor_time, status] = mdsvalue('dim_of(\it)');
    if (mod(status, 2) == 1);
      itf = mdsvalue('\it'); % in amps
      btor = (4*pi*1e-7) * itf *(16*130) /(2*pi*1.8); % about 4327 amps/tesla
      btor = interp1(btor_time, btor, timebase_column);
    else;
      btor = NaN(length(timebase), 1);
    end;
    mdsclose;
  else;
    btor = NaN(length(timebase), 1);
  end;
elseif (shot >= 60000);
  [~, status] = mdsopen('pcs_east', double(shot));
  if (mod(status,2) == 1);
    [btor_time, status] = mdsvalue('dim_of(\it)');
    if (mod(status, 2) == 1);
      itf = mdsvalue('\it'); % in amps
      btor = (4*pi*1e-7) * itf *(16*130) /(2*pi*1.8); % about 4327 amps/tesla
      btor = interp1(btor_time, btor, timebase_column);
    else;
      btor = NaN(length(timebase), 1);
    end;
    mdsclose;
  else;
    btor = NaN(length(timebase), 1);
  end;
else;   % get itf from eng_tree for shots < 60000
  [~, status] = mdsopen('eng_tree', double(shot));
  if (mod(status, 2) == 1);
    [itf, status] = mdsvalue('\it'); % in amps
    if (mod(status, 2) == 1);
      itf = mean(itf);  % Btor is constant in time (superconducting magnet)
      btor = (4*pi*1e-7) * itf *(16*130) /(2*pi*1.8); % about 4327 amps/tesla
      btor = [btor; btor];   % construct 2-point signal from scalar value
      btor_time = [0; 1000]; % construct 2-point timebase
      btor = interp1(btor_time, btor, timebase_column);
    else;
      btor = NaN(length(timebase), 1);
    end;
    mdsclose;
  else;
    btor = NaN(length(timebase), 1);
  end;
end;

% Now get the Mirnov signal(s)

[~, status] = mdsopen('east', shot);
if (mod(status,2) == 0);
  return;
end;

time_window = 0.001;

[time, status] = mdsvalue(['dim_of(\' Mirnov_sensor ')']);
bp_dot = mdsvalue(['\' Mirnov_sensor]);  % tesla/s
mdsclose;
if mod(status,2) == 0 || ...
   length(bp_dot) == 1 || length(time) == 1 || ...
   length(find(~isnan(time)))   <= 1 || ...
   length(find(~isnan(bp_dot))) <= 1 || ...
   length(time) ~= length(bp_dot);
  return; 
end;
for i=1:ntimes
  indices = find(time > (timebase(i)-time_window) & time < timebase(i));
  Mirnov_std(i) = std(bp_dot(indices),'omitnan');
end;
Mirnov_std = transpose(Mirnov_std);

Mirnov_std_normalized = Mirnov_std ./ abs(btor);

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Mirnov_std = transpose(Mirnov_std);
  Mirnov_std_normalized = transpose(Mirnov_std_normalized);
end;

end
