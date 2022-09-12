%retrieve_all_params;
define_indices;
get_ip_sign;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');

delta = 5; % 5 kA bin size
bins = [-250 : delta : 250];

subset = intersect(indices_flattop_disrupt_in_flattop, ...
  find((time_until_disrupt > 0.0095) & (time_until_disrupt < 0.0155)));

slice_data = ip_error(subset)/1.e3 .* ipsign(subset);

counts = hist(slice_data, bins);
bar_handle = bar(bins, counts/sum(counts)/delta, 'hist');
%set(bar_handle, 'EdgeColor', 'none');
%set(bar_handle, 'FaceColor', 'k');
set(gca, 'linewidth', 1.5, 'fontsize', 12)

xlim([-200, 200]);
xlabel('Ip error [kA]', 'FontSize', 14);
ylabel('Probability [kA^{-1}]', 'FontSize', 14);
title('C-Mod: Distribution of Ip error', 'FontSize', 16);

text(0.20, 0.89, '10-15 ms', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.82, 'before', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.75, 'disruption', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.68, '(flattop only)', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('ip_error_disrupt_030_histogram.png', '-dpng');

ylim([0,.035]);
print('ip_error_disrupt_030_histogram_scaled.png', '-dpng');
