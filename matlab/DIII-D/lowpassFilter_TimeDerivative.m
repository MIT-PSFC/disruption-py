function [lowpass_ip,dip_lowpass] = lowpassFilter_TimeDerivative(ip,iptime);
% design a low-pass filter at 50Hz

N   = 100;        % FIR filter order
Fp  = 50;         % 50 Hz passband-edge frequency
Fs  = 2e3;        % 2 kHz sampling frequency
Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
Rst = 1e-4;       % Corresponds to 80 dB stopband attenuation

eqnum = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge'); % eqnum = vec of coeffs
%fvtool(eqnum,'Fs',Fs,'Color','White') % Visualize filter

lowpassFilt =dsp.LowpassFilter('DesignForMinimumOrder',false, ...
    'FilterOrder',N,'PassbandFrequency',Fp,'SampleRate',Fs,...
    'PassbandRipple',0.01, 'StopbandAttenuation',80);

lowpass_ip = lowpassFilt(ip);
dip_lowpass = gradient(lowpass_ip,iptime);
figure()
plot(iptime,dip_lowpass)
xlabel('time [s]')
ylabel('dI_p/dt (lowpass) [A/s]')
end