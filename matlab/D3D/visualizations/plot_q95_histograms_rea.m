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

bins = [0 : .1 : 10.5];
histogram( q95(indices_class_no_disrupt), bins, ...
  'normalization', 'probability', 'displaystyle', 'stairs', ...
  'edgecolor', 'b', 'linewidth', 2);
xlim([1,9]);
ylim([0,.3]);
set(gca, 'fontsize', 18);
xlabel('q_{95}', 'fontsize', 20);
ylabel('Probability', 'fontsize', 20);
%title('Histograms of q_{95}', 'fontsize', 16);

% Overlay histogram of q95 for the disrupt class

hold on
histogram( q95(indices_class_disrupt), bins, ...
  'normalization', 'probability', 'displaystyle', 'stairs', ...
  'edgecolor', 'r', 'linewidth', 2);
%legend('non-disruptions, flattop data', ...
 % 'flattop disruptions (C. Rea''s list), flattop data');
hold off
