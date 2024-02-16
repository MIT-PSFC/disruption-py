mdsconnect('atlas.gat.com');
shot =160942;
shot=161396;
[shotopened, status]=mdsopen('d3d', shot);
ip = mdsvalue(['ptdata("ip", ' num2str(shot) ')']);
[iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' num2str(shot) '))']);
    
if (mod(status,2) == 1);
    iptime = iptime/1.e3; % convert ms to s
    dip_dt = gradient(ip, iptime);
end

[V_loop, status] = mdsvalue(['ptdata("vloopb", ' num2str(shot) ')']); % Volts
[V_loop_time, status] = ...
   mdsvalue(['dim_of(ptdata("vloopb", ' num2str(shot) '))']); % ms
V_loop_time = V_loop_time/1000; % convert to seconds

V_loop_med11 = medfilt1(V_loop, 11); % median filter;  Cristina picked the value '20'
V_loop_med21 = medfilt1(V_loop, 21);
V_loop_pau = my_medfilt1(V_loop,11);
%if (mod(status,2) == 0 || length(V_loop_time) <= 1);
%  V_loop_med = NaN(size(timebase));
%  return;
%else;
%  V_loop_med = interp1(V_loop_time, V_loop_med, iptime, 'linear');
%end;

figure()
plot(V_loop_time,V_loop,'.-b')
hold on
plot(V_loop_time,V_loop_med11,'.-g')
hold on
plot(V_loop_time,V_loop_pau,'.-r')
hold off
xlim([0,8])
%xlim([5.97,5.979])

%calculate disrtortion as rms deviation
dist1 = sqrt((1/length(V_loop_med11)*sum((V_loop-V_loop_med21).^2)));
% coeff. of determination R
%https://en.wikipedia.org/wiki/Coefficient_of_determination
SSres=sum((V_loop-V_loop_pau).^2);
SStot = var(V_loop_pau)*length(V_loop_pau);
R2=1-SSres/SStot;

ip_med11 = medfilt1(ip,11);
ip_med21 = medfilt1(ip,21);
figure()
plot(iptime,ip,'.-b')
hold on
plot(iptime,ip_med11,'.-g')
hold on
plot(iptime,ip_med21,'.-r')
hold off
xlim([0,8])

gsastd_ip10 = GSASTD(iptime,ip_med11,1,16,3,1,0);
gsastd_ip15 = GSASTD(iptime,ip_med11,1,20,3,1,0);
gsastd_ip = GSASTD(iptime,ip,1,20,3,1,0);

figure()
%plot(iptime,gradient(ip,iptime),'.-k')
%hold on
%plot(iptime,gsastd_ip10,'.-b')
plot(iptime,gsastd_ip,'.-b')
hold on
plot(iptime,gsastd_ip15,'.-r')
hold off
xlim([0,8])


%plot(iptime,GSASTD(iptime,median_ip,1,15,3,1,0),'-b')
%hold on
%plot(iptime,ip,'-g')
%hold on
%plot(xlim, [1 1]*(-1.e6),'-r')
%hold off
%xlim([-0.1,8])
