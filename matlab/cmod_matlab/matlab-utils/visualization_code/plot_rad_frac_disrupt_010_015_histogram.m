%retrieve_all_params;
define_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');
delta = 0.01; % bin size
bins = [-1 : delta : 3];

subset = intersect(indices_flattop_disrupt_in_flattop, ...
  find((time_until_disrupt > 0.0095) & (time_until_disrupt < 0.0155)));

slice_data = radiated_fraction(subset);

counts = hist(slice_data, bins);
bar_handle = bar(bins, counts/sum(counts)/delta, 'hist');
%set(bar_handle, 'EdgeColor', 'none');
%set(bar_handle, 'FaceColor', 'k');
set(gca, 'linewidth', 1.5, 'fontsize', 12)

xlim([0, 2]);
xlabel('P_{rad}/P_{input}', 'FontSize', 14);
ylabel('Histogram', 'FontSize', 14);
title('C-Mod: Histogram of P_{rad}/P_{input}', 'FontSize', 16);

text(0.80, 0.89, '10-15 ms', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.82, 'before', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.75, 'disruption', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.68, '(flattop only)', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('east_rad_frac_disrupt_010-015_histogram.png', '-dpng');

ylim([0,3]);
print('east_rad_frac_disrupt_010-015_histogram_scaled.png', '-dpng');
