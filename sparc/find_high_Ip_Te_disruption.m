%% Retrieve from database
%{
db = set_database('logbook');
result = fetch(db,['select shot, time, time_until_disrupt, dipprog_dt, ip, ' ...
		'beta_n, n_e, ip_error from disruption_warning order by shot, time']);

shot = cell2mat(result(:,1));
time = cell2mat(result(:,2));
time_until_disrupt = cell2mat(result(:,3));
dipprog_dt = cell2mat(result(:,4));
ip = cell2mat(result(:,5));
beta_n = cell2mat(result(:,6));
n_e = cell2mat(result(:,7));
ip_error = cell2mat(result(:,8));
ip_prog = ip-ip_error;

define_indices

%% Plot distributions of ip, beta_n/n_e
figure()
histogram(beta_n./n_e,'normalization','probability')

%% Plot distributions of ip, beta_n/n_e
figure()
histogram(abs(ip),'normalization','probability')

%% Retrieve list of shots with high beta_n/ne (> 150) and high Ip (> 1e6)
shotlist = unique(cell2mat(fetch(db,['select shot from disruption_warning where ' ...
		'ip>1e6 and time_until_disrupt is not null and beta_n > 1.2'])));
%}
% Plot individual shot
indices = find(shot==1160922032);
figure(); hold on;
plot(time(indices),ip(indices)/1e6,'k-','DisplayName','$I_p$')
plot(time(indices),ip_prog(indices)/1e6,'k--','DisplayName','$I_{p,prog}$')
xlabel('Time [s]')
ylabel('Current [MA]')
lgd = legend();
lgd.('Interpreter') = 'latex';
