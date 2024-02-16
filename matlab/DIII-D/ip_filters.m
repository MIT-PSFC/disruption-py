mdsconnect('atlas.gat.com');
shot =160942;
WindowLength = 8;
    
[shotopened, status]=mdsopen('d3d', shot);
ip = mdsvalue(['ptdata("ip", ' num2str(shot) ')']);
[iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' num2str(shot) '))']);
    
if (mod(status,2) == 1);
    iptime = iptime/1.e3; % convert ms to s
    dip_dt = gradient(ip, iptime);
end

medFilt = dsp.MedianFilter(WindowLength);
movavgWin = dsp.MovingAverage(WindowLength);
ip_mvavg = movavgWin(ip);
ip_med = medFilt(ip);
    
Fs = 1./(iptime(2)-iptime(1));%length(iptime);% sample rate in Hertz
    
  %{
    scope  = dsp.TimeScope('SampleRate',Fs,...
    'TimeSpanOverrunAction','Scroll',...
    'TimeDisplayOffset',iptime,'ShowGrid',true,...%'TimeSpan',1,'ShowGrid',true,...
    'TimeSpanSource', 'Auto',...
    'YLimits',[-2,12]*10^5,...
    'LayoutDimensions',[3 1],...
    'NumInputPorts',3);

    scope.ActiveDisplay = 1;
    scope.Title = 'Ip [MA]';
    scope.ActiveDisplay = 2;
    scope.Title = 'Moving Average Output (Window Length = 10)';
    scope.ActiveDisplay = 3;
    scope.Title = 'Median Filter Output (Window Length = 10)';
    
    
    scope(ip,ip_mvavg,ip_med);
%}
scope2=dsp.TimeScope('SampleRate',Fs,...
       'TimeSpanOverrunAction','Scroll',...
       'TimeDisplayOffset',iptime,'ShowGrid',true,...%'TimeSpan',1,'ShowGrid',true,...
       'TimeSpanSource', 'Auto',...
       'YLimits',[-5,5]*10^6,...
       'LayoutDimensions',[4 1],...
       'NumInputPorts',4);
scope2.ActiveDisplay = 1;
scope2.Title = 'dIp_dt [MA/s]';
scope2.ActiveDisplay = 2;
scope2.Title = 'Median Filter Output 10-th order';
scope2.ActiveDisplay = 3;
scope2.Title = 'GSASTD Output (Window Length = )';
scope2.ActiveDisplay = 4;
scope2.Title = 'Low-pass Time Derivative Output (Window Length = )';

median_ip = medfilt1(ip,10);
dip_median = gradient(median_ip,iptime);
gsastd_ip = GSASTD(iptime,median_ip,1,10,3,1,0);
%[lowpass_ip,dip_lowpass] = lowpassFilter_TimeDerivative(ip,iptime);
scope2(dip_dt, dip_median, gsastd_ip, dip_lowpass);
%end
plot(iptime,GSASTD(iptime,median_ip,1,15,3,1,0),'-b')
hold on
plot(iptime,ip,'-g')
hold on
plot(xlim, [1 1]*(-1.e6),'-r')
hold off
xlim([-0.1,8])
