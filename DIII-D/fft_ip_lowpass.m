mdsconnect('atlas.gat.com');
shot =160942;
[shotopened, status]=mdsopen('d3d', shot);
ip = mdsvalue(['ptdata("ip", ' num2str(shot) ')']);
[iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' num2str(shot) '))']);
 
if (mod(status,2) == 1);
    iptime = iptime/1.e3; % convert ms to s
    %dip_dt = gradient(ip, iptime);
end

Fs = 1/(iptime(2)-iptime(1)); %sampling freq: 2kHz
T = 1/Fs;
L = length(ip);
% fft of ip signal
ip_fft = fft(ip);
P2 = abs(ip_fft/L);
P1 = P2(1:L/2+1);
%excluding the zero-frequency (DC) component
P1(2:end-1) = 2*P1(2:end-1);
%define frequency domain and plot the single-sided spectrum
f = Fs*(0:(L/2))/L;
plot(f,P1)
title('single sided ampl spectrum of Ip(t), shot 160942')
xlabel('f [Hz]')
ylabel('|P1(f)|')
%xlim([-0.1,10])
ylim([0,1000])

% design a low-pass filter at 50Hz

N   = 100;        % FIR filter order
Fp  = 50;         % 50 Hz passband-edge frequency
Fs  = 2e3;        % 2 kHz sampling frequency
Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
Rst = 1e-4;       % Corresponds to 80 dB stopband attenuation

eqnum = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge'); % eqnum = vec of coeffs
fvtool(eqnum,'Fs',Fs,'Color','White') % Visualize filter

lowpassFilt =dsp.LowpassFilter('DesignForMinimumOrder',false, ...
    'FilterOrder',N,'PassbandFrequency',Fp,'SampleRate',Fs,...
    'PassbandRipple',0.01, 'StopbandAttenuation',80);

lowpass_ip = lowpassFilt(ip);
dip_lowpass = gradient(lowpass_ip,iptime);
figure()
plot(iptime,dip_lowpass)
xlabel('time [s]')
ylabel('dI_p/dt (lowpass) [A/s]')