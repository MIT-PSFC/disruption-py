db = set_database_d3d('d3drdb');
retrieve_all_params_;
define_indices;

indices_2015 = find(shot>1150000000 & shot<1160000000);

ii = intersect(indices_flattop_no_disrupt,indices_2015);
jj = intersect(indices_flattop_disrupt_in_flattop,indices_2015);
jj_far = intersect(jj, find(time_until_disrupt > 0.040));
jj_near = intersect(jj, find(time_until_disrupt <= 0.040 & ...
  time_until_disrupt > 0.002)); % 2 ms 'black window'

bins = 0.00 : 0.02 : 2.00;

histogram(li(ii), bins, 'normalization', 'probability', ...
  'facecolor', 'blue');
xlim([0.5, 2.0]);
ylim([0.00, 0.11]);
set(gca, 'fontsize', 12);
xlabel('l_i', 'fontsize', 14);
ylabel('probability histogram', 'fontsize', 14);
title('Histograms of l_i for Alcator C-Mod', 'fontsize', 16);
hold on;
histogram(li(jj_far), bins, 'normalization', 'probability', ...
  'facecolor', 'green');
histogram(li(jj_near), bins, 'normalization', 'probability', ...
  'facecolor', 'red');
legend('non-disruptive', 'time until disrupt > 40 ms', ...
  'time until disrupt < 40 ms', 'location', 'northwest');
hold off;