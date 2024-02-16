function [n1rms, n2rms, n1rms_normalized, n2rms_normalized] = get_n1rms_n2rms(shot, timebase)

% Read in the saddle sensor data and the rmp currents from the MDSplus
% tree.  All the outputs have time as their 1st dimension,
% i.e. rmp(time,rmpcoil#) and sad(time,saddlecoil#), and the time arrays
% are column vectors.

% The ordering of the rmp data is: irmpu1 to irmpu8, then irmpl1 to irmpl8.
% The ordering of the saddle data is: sad_pa, sad_bc, .. sad_no.
if (size(timebase,2) > 1);
  timebase = transpose(timebase);
end;


ntimes = length(timebase);
n1rms = NaN(size(timebase));
n1rms_normalized = NaN(size(timebase));
n2rms = NaN(size(timebase));
n2rms_normalized = NaN(size(timebase));


mdsconnect('mds.ipp.ac.cn');

[~, status] = mdsopen('east', shot);
if mod(status,2) == 0;
  mdsclose;
  return;
end;




[Mirtime, status] = mdsvalue('dim_of(\mitab2)');
if mod(status,2) == 1 && length(Mirtime) > 2;
  Mir = NaN(length(Mirtime), 16);

  Mir(:,1) = mdsvalue('\mitab2');
  Mir(:,2) = mdsvalue('\mitbc2');
  Mir(:,3) = mdsvalue('\mitcd2');
  Mir(:,4) = mdsvalue('\mitde2');
  Mir(:,5) = mdsvalue('\mitef2');
  Mir(:,6) = mdsvalue('\mitfg2');
  Mir(:,7) = mdsvalue('\mitgh2');
  Mir(:,8) = mdsvalue('\mithi2'); 
  Mir(:,9) = mdsvalue('\mitij2');
  Mir(:,10) = mdsvalue('\mitjk2');
  Mir(:,11) = mdsvalue('\mitkl2');
  Mir(:,12) = mdsvalue('\mitlm2');
  Mir(:,13) = mdsvalue('\mitmn2');
  Mir(:,14) = mdsvalue('\mitno2');
  Mir(:,15) = mdsvalue('\mitop2');
  Mir(:,16) = mdsvalue('\mitpa2');
else
    mdsclose;
    return
end
mdsclose;

if (shot > 65165);
  [~, status] = mdsopen('eng_tree', double(shot));
  if (mod(status, 2) == 1);
    [btor_time, status] = mdsvalue('dim_of(\it)');
    if (mod(status, 2) == 1);
      itf = mdsvalue('\it'); % in amps
      btor = (4*pi*1e-7) * itf *(16*130) /(2*pi*1.8); % about 4327 amps/tesla
      btor = interp1(btor_time, btor, timebase);
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
      btor = interp1(btor_time, btor, timebase);
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
      btor = interp1(btor_time, btor, timebase);
    else;
      btor = NaN(length(timebase), 1);
    end;
    mdsclose;
  else;
    btor = NaN(length(timebase), 1);
  end;
end;

output = fft(Mir, [], 2);
amplitude = abs(output)/size(output,2);
amplitude(:, 2:end) = amplitude(:, 2:end) * 2;
phase = atan2(imag(output), real(output));
n1=amplitude(:,2).*cos(phase(:,2));
n2=amplitude(:,3).*cos(phase(:,3));

time_window = 0.001;

for i=1:ntimes
  indices = Mirtime > (timebase(i)-time_window) & Mirtime < (timebase(i)+time_window);
  n1rms(i) = std(n1(indices),'omitnan');
  n2rms(i) = std(n2(indices),'omitnan');
end;


n1rms_normalized = n1rms ./ abs(btor);
n2rms_normalized = n2rms ./ abs(btor);


end
