retrieve_all_params;
define_indices;

z_norm = zcur(indices_flattop_disrupt_in_flattop); 
t = time_until_disrupt(indices_flattop_disrupt_in_flattop);
idx = find(abs(z_norm)>4.5*1e-2 & t<0.1);

% idx = find(abs(z_norm)>2*1e-2 & t<0.1 & t>0.03); gives 35 shots

s = shot(indices_flattop_disrupt_in_flattop);
shotlist = unique(s(idx));

100*length(shotlist)/length(unique(s))%17.4% = 217/1244

%idx = find(abs(z_norm)>4.5*1e-2 & t>=0.1 & t<0.350); just 157790

idx = find(abs(z_norm)<=2*1e-2 & t<=0.03); % anything stable during CQ? 1240 shots

ii = find(shot==176370);
plot(time(ii),ip(ii)*1e-6) 
xlabel('time [s]')
ylabel('Ip [MA]')
title('176370')
xlim([3.3335,3.355])

figure()
plot(time(ii),zcur(ii)*1e-2)
xlabel('time [s]')
ylabel('Z\_cur [cm]')
title('176370')
xlim([3.3335,3.355])


ii = find(shot==175353);%error field locking
plot(time(ii),ip(ii)*1e-6)
xlabel('time [s]')
ylabel('Ip [MA]')
title('175353')
xlim([3.34,4.47])

figure()
plot(time(ii),zcur(ii)*1e2)
xlabel('time [s]')
ylabel('Z\_cur [cm]')
title('175353')
xlim([3.34,4.47])

