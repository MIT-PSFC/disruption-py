function [rmptime, rmp, sadtime, sad] = get_rmp_and_saddle_signals(shot);

% Read in the saddle sensor data and the rmp currents from the MDSplus
% tree.  All the outputs have time as their 1st dimension,
% i.e. rmp(time,rmpcoil#) and sad(time,saddlecoil#), and the time arrays
% are column vectors.

% The ordering of the rmp data is: irmpu1 to irmpu8, then irmpl1 to irmpl8.
% The ordering of the saddle data is: sad_pa, sad_bc, .. sad_no.

mdsconnect('mds.ipp.ac.cn');
[shotopened, status] = mdsopen('east', shot);
if mod(status,2) == 0;
  rmptime = [];
  rmp = [];
  sadtime = [];
  sad = [];
  return;
end;

[rmptime, status] = mdsvalue('dim_of(\irmpu1)');
if mod(status,2) == 1 && length(rmptime) > 2;
  rmp = NaN(length(rmptime), 16);

  signal = mdsvalue('\irmpu1');
  if (length(signal) == length(rmptime));
    rmp(:, 1) = signal;
  else;
    rmp(:, 1) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpu2');
  if (length(signal) == length(rmptime));
    rmp(:, 2) = signal;
  else;
    rmp(:, 3) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpu3');
  if (length(signal) == length(rmptime));
    rmp(:, 3) = signal;
  else;
    rmp(:, 3) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpu4');
  if (length(signal) == length(rmptime));
    rmp(:, 4) = signal;
  else;
    rmp(:, 4) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpu5');
  if (length(signal) == length(rmptime));
    rmp(:, 5) = signal;
  else;
    rmp(:, 5) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpu6');
  if (length(signal) == length(rmptime));
    rmp(:, 6) = signal;
  else;
    rmp(:, 6) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpu7');
  if (length(signal) == length(rmptime));
    rmp(:, 7) = signal;
  else;
    rmp(:, 7) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpu8');
  if (length(signal) == length(rmptime));
    rmp(:, 8) = signal;
  else;
    rmp(:, 8) = NaN(size(rmptime));
  end;

  signal = mdsvalue('\irmpl1');
  if (length(signal) == length(rmptime));
    rmp(:, 9) = signal;
  else;
    rmp(:, 9) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpl2');
  if (length(signal) == length(rmptime));
    rmp(:,10) = signal;
  else;
    rmp(:,10) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpl3');
  if (length(signal) == length(rmptime));
    rmp(:,11) = signal;
  else;
    rmp(:,11) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpl4');
  if (length(signal) == length(rmptime));
    rmp(:,12) = signal;
  else;
    rmp(:,12) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpl5');
  if (length(signal) == length(rmptime));
    rmp(:,13) = signal;
  else;
    rmp(:,13) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpl6');
  if (length(signal) == length(rmptime));
    rmp(:,14) = signal;
  else;
    rmp(:,14) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpl7');
  if (length(signal) == length(rmptime));
    rmp(:,15) = signal;
  else;
    rmp(:,15) = NaN(size(rmptime));
  end;
  signal = mdsvalue('\irmpl8');
  if (length(signal) == length(rmptime));
    rmp(:,16) = signal;
  else;
    rmp(:,16) = NaN(size(rmptime));
  end;

%{
  rmp(:, 1) = mdsvalue('\irmpu1');
  rmp(:, 2) = mdsvalue('\irmpu2');
  rmp(:, 3) = mdsvalue('\irmpu3');
  rmp(:, 4) = mdsvalue('\irmpu4');
  rmp(:, 5) = mdsvalue('\irmpu5');
  rmp(:, 6) = mdsvalue('\irmpu6');
  rmp(:, 7) = mdsvalue('\irmpu7');
  rmp(:, 8) = mdsvalue('\irmpu8');

  rmp(:, 9) = mdsvalue('\irmpl1');
  rmp(:,10) = mdsvalue('\irmpl2');
  rmp(:,11) = mdsvalue('\irmpl3');
  rmp(:,12) = mdsvalue('\irmpl4');
  rmp(:,13) = mdsvalue('\irmpl5');
  rmp(:,14) = mdsvalue('\irmpl6');
  rmp(:,15) = mdsvalue('\irmpl7');
  rmp(:,16) = mdsvalue('\irmpl8');
%}
else;
  rmptime = [];
  rmp = [];
end;

[sadtime, status] = mdsvalue('dim_of(\sad_pa)');
if mod(status,2) == 1 && length(sadtime) > 2;
  sad = NaN(length(sadtime), 8);

  sad(:,1) = mdsvalue('\sad_pa');
  sad(:,2) = mdsvalue('\sad_bc');
  sad(:,3) = mdsvalue('\sad_de');
  sad(:,4) = mdsvalue('\sad_fg');
  sad(:,5) = mdsvalue('\sad_hi');
  sad(:,6) = mdsvalue('\sad_jk');
  sad(:,7) = mdsvalue('\sad_lm');
% sad(:,8) = mdsvalue('\sad_no'); % not operational in 2015
  sad(:,8) = mdsvalue('\sad_lo') - mdsvalue('\sad_lm');
else;
  sadtime = [];
  sad = [];
end;

mdsclose;
return;
end
