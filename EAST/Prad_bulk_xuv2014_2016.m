% The original of this program was given to us by Duan Yanmin
% <ymduan@ipp.ac.cn>.  I (R. Granetz) modified it to reduce the smoothing
% window.  Duan Yanmin uses a 100 ms window (window width = 100 points, and
% the data is sampled at 1 kHz).  This is too much filtering when studying
% disruptions.  I also made other modifications, including returning the
% Prad signal and its timebase, as well as an error status parameter, and I
% commented out the 'plot' statement.  June 2016
% ----------------

%%Program to caculate total radiated power in bulk plasma. Unit/kW
%%Just for axuv datas from 2014 -2016 experimental campaign.

% I changed units to watts -- R. Granetz

%function shotno1=Prad_bulk_xuv2014_2016(shotno0)
function [data, time, status] = Prad_bulk_xuv2014_2016(shotno0)
%%%%%%%%%%%%%%%%%%
%  clear all
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%calibration factors%%%
Fac1=[1.3681,1.3429,1.3215,1.3039,1.2898,1.2793,1.2723,1.2689,1.2689,1.2723,1.2793,1.2898,1.3039,1.3215,1.3429,1.3681]*10.^(4);
Fac2=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];% factors of Amp.Gain
Fac3=[1,1,1,1]; %cross calibration factors between arrays
Fac4=1*10.^(-3);%unit convert
Fac5=2.5;% corrected factor by cross calibration with foil bolometer
Maj_R=1.85;
Del_r=[3.6,3.6,3.5,3.4,3.4,3.3,3.3,3.2,3.2,3.1,3.1,3.0,3.0,2.9,2.9,2.8,...
  2.9,2.8,2.8,2.8,2.8,2.8,2.8,2.7,2.7,2.7,2.7,2.7,2.6,2.6,2.6,2.6,...
  2.6,2.6,2.6,2.6,2.7,2.7,2.7,2.7,2.7,2.8,2.8,2.8,2.8,2.8,2.8,2.9,...
  2.8,2.9,2.9,3.0,3.0,3.1,3.1,3.2,3.2,3.3,3.3,3.4,3.4,3.5,3.6,3.6]*0.01;

%%%%%%%%readdata%%%%%%%
%mdsconnect('202.127.204.12');
mdsconnect('mds.ipp.ac.cn');
[shotnumber, status_mdsopen] = ...
  mdsopen('EAST_1',double(shotno0));% tree:east_1,sampling fre.1KHz
if (mod(status_mdsopen, 2)==0);
  status = 0;
  data = NaN;
  time = NaN;
  return;
end;
%%%%%%%%%%%%%%%%%
for j=1:4
for i=1:16  
    m=(j-1)*16+i;
    istr=num2str(m); 
%------ These modifications by R. Granetz, June 2016 -----
    [data1, status_mdsvalue] = mdsvalue(strcat('\pxuv',istr));
    if (mod(status_mdsvalue, 2) == 0);
      status = 0;
      data = NaN;
      time = NaN;
      return;
    end;
%---------------------------------------------------------
    data1=data1(:);
    num=size(data1);
    nn=num(1);
    xuv(1:nn,m)=data1; % Note: the array "xuv" is not pre-dimensioned?
 xuv(1:nn,m)=xuv(1:nn,m)-mean(xuv(1:100,m)); % substract Zero drift 

%------ this modification by R. Granetz, June 2016 -----
%xuv(1:nn,m)=smooth(xuv(1:nn,m),100);  %%smooth function
 timebase=mdsvalue('dim_of(\pxuv32)');
 dt = timebase(2)-timebase(1); % sampling interval in seconds
 smoothing_time = 1.e-3; % seconds
 smoothing_window = max([round(smoothing_time/dt), 1]);
 xuv(1:nn,m)=smooth(xuv(1:nn,m),smoothing_window);  %%smooth function
%-------------------------------------------------------

 xuv(1:nn,m)=xuv(1:nn,m)*Fac1(i)*Fac2(j,i)*Fac4*Fac3(j);
end
end
%%%%%correction for bad channels
 xuv(:,12)=0.5*(xuv(:,11)+xuv(:,13));
 xuv(:,36)=0.5*(xuv(:,35)+xuv(:,36));
 %%%%%%
time = mdsvalue('dim_of(\pxuv32)'); % [s]
mdsclose;
%%%%%%%%%%%%
temp1=zeros(nn,64);
for i=1:64
    temp1(:,i)=xuv(:,i)*2*pi*Maj_R*Del_r(i);
end
 temp2=[temp1(:,8:15) temp1(:,21:32) temp1(:,33:48) temp1(:,54:59)];% substract divertor measuring and overlap measuring
 Prm=sum(temp2,2)*Fac5;
%%%%%%%%%
data = Prm;    %%unit/kW, the total radiated power in bulk plasma
data = data * 1.e3; % convert to watts (R. Granetz)
status = 1;
%figure(100)
%plot(time,data);
