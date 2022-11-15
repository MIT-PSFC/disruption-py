function a = getbolo(ishot, drtau);

if ~exist('drtau', 'var'); drtau = 50; end;
drtau = drtau/1.e3;
gam=zeros(1,49);
tau=zeros(1,49);
kappa=[1.2307e8,1.4719e8,1.5750e8,1.5316e8,1.6025e8,1.6418e8, ...
  1.7073e8,2.3007e8,2.2421e8,2.1928e8,2.1501e8,2.1228e8,2.0982e8, ...
  2.0762e8,2.0568e8,1.8170e8,1.4038e8,1.1152e8,0.98347e8,0.91370e8, ...
  0.86452e8,0.87789e8,0.93361e8,0.97920e8,5.202e8,5.072e8,4.865e8, ...
  4.527e8,4.180e8,4.154e8,4.137e8,4.248e8,4.459e8,4.727e8,5.195e8, ...
  6.633e7,6.544e7,6.699e7,7.100e7,1.749e8,1.860e8,1.939e8,2.016e8, ...
  2.095e8,2.169e8,2.256e8,2.311e8,10.39e7,4.4559e8,4.3989e8,4.1525e8, ...
  3.9185e8,3.6744e8,3.4002e8,3.1598e8];
scrfact=ones(1,48);

if ishot > 97000;
 kappa(25:35)=[1.5398e8,1.5013e8,1.4461e8,1.3399e8,1.2371e8,1.2297e8, ...
  1.2246e8,1.2574e8,1.3198e8,1.3993e8,1.5376e8];
% transmission factor through ECH screen
%
%improved ECH screen
 scrfact=[0.7990,0.8006,0.8021,0.8035,0.8048,0.8060, ...
 	  0.8070,0.8080,0.8089,0.8099,0.8105,0.8096, ...
	  0.8080,0.8058,0.8035,0.8011,0.8038,0.8062, ...
	  0.8083,0.8100,0.8096,0.8079,0.8061,0.8042, ...
	  0.8051,0.8062,0.8074,0.8085,0.8096,0.8101, ...
	  0.8079,0.8056,0.8031,0.8004,0.7973,0.7991, ...
	  0.8019,0.8044,0.8066,0.8086,0.8096,0.8105, ...
	  0.8099,0.8089,0.8080,0.8070,0.8060,0.8037];
end;

onechannel=struct('label', '', ...
                  'R', 0.0, ...
                  'Z', 0.0, ...
                  'ier', int32(0), ...
                  'pwr', zeros(1,4096), ...
                  'raw', zeros(1,16384), ...
                  'gam', 0.0, ...
                  'tau', 0.0, ...
                  'scrfact', 0.0);
ch=repmat(onechannel,1,48);

a=struct('shot', int32(0), ...
         'kappa', kappa, ...
         'time', zeros(1,4096), ...
         'rawtime', zeros(1,16384), ...
         'ntimes', int32(0), ...
         'totpwr', zeros(1,4096), ...
         'chan', ch);

mdsconnect('atlas.gat.com');
[~, status] = mdsopen('bolom', ishot);

ier = int32(0);

chname = cell(48,1);
for i=1:24;
  chname{i}=['bol_u' num2str(i, '%02i') '_v'];
  chname{i+24}=['bol_l' num2str(i, '%02i') '_v'];
end;

a.shot = ishot;
if ishot > 79400;
  data = mdsvalue('\bol_prm');
  gam = data(1:49);
  tau = data(50:98);
else;
end;

time = mdsvalue(['dim_of(\top.raw:' chname{1} ')']) / 1.e3;

% The bolometry signals and their common timebase should each have 163840
% samples, digitized at 10 kHz.  Tony Leonard's IDL procedure down-samples
% these to 16384 samples using IDL's REBIN function, and then again to just
% 4096 samples.  But REBIN does a forward averaging, which introduces some
% amount of non-causality.  Instead, for the timebase I will generate a new
% array with 16384 points, covering the same range of time.  For the
% bolometry signals, I will use nearest-neighbor linear interpolation to
% down-sample.

% Sometimes the bolo data is garbage.  I need to execute a few checks to
% decide if the data is valid or not.  If not, set 'ier' = 1 in the
% appropriate field for each channel.

if (length(time) <= 16384 || ...
    time(end) <= time(1)  || ...
    min(diff(time)) < -1e-5);
  for i = 1:48;
    a.chan(i).ier = int32(1);
  end;
  return;
end;

time = linspace(min(time), max(time), 16384);
dt = time(2)-time(1);
windowsize = round(drtau / dt);
smoothing_kernal = (1/windowsize) * ones(1, windowsize);

a.ntimes = int32(length(time)/4);
a.time = linspace(min(time), max(time), a.ntimes);
tdel = a.time(2)-a.time(1);
a.rawtime = time;
m = 2*fix(fix(1000*drtau)/fix(1000*tdel)/2)+1;
k = [0:m-1] - fix((m-1)/2);
nzer = find(k ~= 0);
k(nzer) = 1.0 ./ k(nzer);
k = k/tdel/(fix(m/2)*2);

for i=1:48;
  a.chan(i).label=chname{i};
  data_orig = mdsvalue(['\top.raw:' chname{i}]);
  time = mdsvalue(['dim_of(\top.raw:' chname{i} ')'])/1.e3;

% Reduce to 16384 points using Matlab's INTERP1 routine with
% nearest-neighbor linear interpolation

  data = interp1(time, data_orig, a.rawtime);

  a.chan(i).ier = int32(0);
  a.chan(i).raw = data;
  a.chan(i).gam = gam(i+1);
  a.chan(i).tau = tau(i+1);
  a.chan(i).scrfact = scrfact(i);

% Subtract baseline offset
  temp = data - sum(data(1:20))/20;

% Filter signal using causal moving average filter (i.e. boxcar)
  temp_filtered = filter(smoothing_kernal, 1, temp);

% Calculate derivative of signal
  drdt = gradient(temp_filtered, dt);

  a.chan(i).pwr = ...
    medfilt1((gam(i+1)*temp_filtered + tau(i+1)*drdt)/scrfact(i), windowsize);
end;

mdsclose;
return;
end
