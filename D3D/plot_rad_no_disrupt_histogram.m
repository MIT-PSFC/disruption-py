%retrieve_all_params;
define_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');
delta = 0.25; % bin size [MW]
bins = [-1 : delta : 16];

subset = indices_flattop_no_disrupt;

slice_data = p_rad(subset)/1.e6;

counts = hist(slice_data, bins);
bar_handle = bar(bins, counts/sum(counts)/delta, 'hist');
%set(bar_handle, 'EdgeColor', 'none');
%set(bar_handle, 'FaceColor', 'k');
set(gca, 'linewidth', 1.5, 'fontsize', 12)

xlim([0, 15]);
xlabel('P_{rad} [MW]', 'FontSize', 14);
ylabel('Histogram [MW^{-1}]', 'FontSize', 14);
title('DIII-D: Histogram of P_{rad}', 'FontSize', 16);

text(0.80, 0.89, 'All', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.82, 'non-disruptive', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.75, 'shots', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.68, '(flattop only)', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('d3d_rad_no_disrupt_histogram.png', '-dpng');


figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');
delta = 0.05; % bin size [MW]
bins = [-1 : delta : 16];

subset = indices_flattop_no_disrupt;

slice_data = radiated_fraction(subset);

counts = hist(slice_data, bins);
bar_handle = bar(bins, counts/sum(counts)/delta, 'hist');
%set(bar_handle, 'EdgeColor', 'none');
%set(bar_handle, 'FaceColor', 'k');
set(gca, 'linewidth', 1.5, 'fontsize', 12)

xlim([0, 2]);
xlabel('P_{rad}/P_{input}', 'FontSize', 14);
ylabel('Histogram', 'FontSize', 14);
title('DIII-D: Histogram of P_{rad}/P_{input}', 'FontSize', 16);

text(0.80, 0.89, 'All', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.82, 'non-disruptive', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.75, 'shots', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.68, '(flattop only)', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('d3d_rad_frac_no_disrupt_histogram.png', '-dpng');
