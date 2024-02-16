db = set_database('d3drdb');
result = fetch(db, ['select shot,time,time_until_disrupt, dipprog_dt, ip,' ...
  'ip_error, power_supply_railed, intentional_disruption, other_hardware_failure,' ...
  'kappa_area,kappa from disruption_warning order by shot,time']);

shot                   = cell2mat(result(:,1));
time                   = cell2mat(result(:,2));
time_until_disrupt     = cell2mat(result(:,3));
dipprog_dt             = cell2mat(result(:,4));
ip                     = cell2mat(result(:,5));
ip_error               = cell2mat(result(:,6));
power_supply_railed    = cell2mat(result(:,7));
intentional_disruption = cell2mat(result(:,8));
other_hardware_failure = cell2mat(result(:,9));
kappa_area             = cell2mat(result(:,10));
kappa                  = cell2mat(result(:,11));

define_indices;

ii = indices_flattop_no_disrupt;
jj = indices_flattop_disrupt_in_flattop;
jj_far = intersect(jj, find(time_until_disrupt > 0.35));
jj_near = intersect(jj, find(time_until_disrupt <= 0.35 & ...
  time_until_disrupt > 0.01)); % 10 ms 'black window'

% 
bins = 1 : 0.01 : 2.1;
histogram(kappa_area(ii), bins, 'normalization', 'probability', ...
  'DisplayStyle', 'stairs', 'EdgeColor', 'g');                                                           
xlim([1,2.1]);
ylabel('Fraction of flattop time slices');
title('DIII-D');
%ylim([0.0,0.08]);
set(gca, 'fontsize', 12);
xlabel('\kappa_{area}', 'fontsize', 14);
hold on;
histogram(kappa_area(jj_far), bins, 'normalization', 'probability', ...
  'DisplayStyle', 'stairs', 'EdgeColor', 'b');
histogram(kappa_area(jj_near), bins, 'normalization', 'probability', ...
  'DisplayStyle', 'stairs', 'EdgeColor', 'r');
legend('non-disruptive', 'time until disrupt > 350 ms', ...
  '10 ms < time until disrupt =< 350 ms', 'location', 'northwest');
hold off;

figure();
bins = 1 : 0.01 : 2.05;
histogram(kappa(ii), bins, 'normalization', 'probability', ...
  'DisplayStyle', 'stairs', 'EdgeColor', 'g');
xlim([1,2]);
ylabel('Fraction of flattop time slices');
title('DIII-D');
%ylim([0.0,0.08]);
set(gca, 'fontsize', 12);
xlabel('\kappa_{sep}', 'fontsize', 14);
hold on;
histogram(kappa(jj_far), bins, 'normalization', 'probability', ...
  'DisplayStyle', 'stairs', 'EdgeColor', 'b');
histogram(kappa(jj_near), bins, 'normalization', 'probability', ...
  'DisplayStyle', 'stairs', 'EdgeColor', 'r');
legend('non-disruptive', 'time until disrupt > 350 ms', ...
  '10 ms < time until disrupt =< 350 ms', 'location', 'northwest');
hold off;
