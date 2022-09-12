cmap_hist = ametrine(256); %Color-blind friendly colormap in /home/montes/matlab
db = set_database('logbook');
result = fetch(db,['select shot, time, time_until_disrupt, dipprog_dt, ip, v_loop, ' ...
	'Greenwald_fraction, beta_p, radiated_fraction from disruption_warning order by shot, time']);
shot = int64(cell2mat(result(:,1)));
time = cell2mat(result(:,2));
time_until_disrupt = cell2mat(result(:,3));
dipprog_dt = cell2mat(result(:,4));
ip = cell2mat(result(:,5)); 
v_loop = cell2mat(result(:,6));
Greenwald_fraction = cell2mat(result(:,7));
beta_p = cell2mat(result(:,8));
radiated_fraction = cell2mat(result(:,9));

% Get sign of current ip
shotlist = unique(shot);
ipsign = nan(size(ip));
for ishot = 1:length(shotlist);
  qq = find(shot == shotlist(ishot));
  plasma_current = ip(qq);
  plasma_current = plasma_current(find(~isnan(plasma_current)));
  ip_integrated = sum(plasma_current);
  ipsign(qq) = ones(size(qq)) * sign(ip_integrated);
end;

% Get indices
define_indices;
ii = indices_flattop_no_disrupt;
jj = indices_flattop_disrupt_in_flattop;
jj_far = intersect(jj,find(time_until_disrupt > 0.35));
jj_near = intersect(jj,find(time_until_disrupt <= 0.04 & time_until_disrupt > 0));

% Generate histogram for parameter of interest
y = v_loop.*ipsign;
bins = 0:0.05:10;

figure(); hold on; colormap(cmap_hist);
histogram(y(ii),bins,'normalization','probability');
histogram(y(jj_far),bins,'normalization','probability');
histogram(y(jj_near),bins,'normalization','probability');
legend({'non-disruptive','far from disruption (> 40 ms)','near disruption (< 40 ms)'},'location','northeast','FontSize',16);
xlabel('$V_{loop}$','Interpreter','latex','FontSize',20);
ylabel('probability','FontSize',20);
xlim([0,5]);
set(gca,'fontsize',18);
