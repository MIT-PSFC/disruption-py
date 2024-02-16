db = set_database('d3drdb');
if ~exist('Greenwald_fraction', 'var');
  retrieve_all_params;
  define_indices;
end;

shotlist=dlmread('../new_flattop_disruptions_cleaned.txt');
shotlist=int32(shotlist);

indices_class_no_disrupt = indices_flattop_no_disrupt;

indx_clean_shots = find(ismember(shot, shotlist));
indices_class_disrupt = ...
  intersect(indices_flattop_disrupt_in_flattop, indx_clean_shots);

% Plot histogram of q95 for the no-disrupt class

bins = [0 : .0001 : 0.005];
histogram( n_equal_1_mode(indices_class_no_disrupt), bins, ...
  'normalization', 'probability', 'displaystyle', 'stairs', ...
  'edgecolor', 'b', 'linewidth', 2);
xlim([0.0001,0.0036]);
ylim([0,0.02]);
set(gca, 'fontsize', 12);
xlabel('n=1 amplitude', 'fontsize', 14);
ylabel('Probability', 'fontsize', 14);
title('Histograms of n=1 amplitude', 'fontsize', 16);

% Overlay histogram of q95 for the disrupt class

hold on
histogram( n_equal_1_mode(indices_class_disrupt), bins, ...
  'normalization', 'probability', 'displaystyle', 'stairs', ...
  'edgecolor', 'r', 'linewidth', 2);
legend('non-disruptions, flattop data', ...
  'flattop disruptions (C. Rea''s list), flattop data');
hold off
