%retrieve_all_params;
define_indices;
get_ip_sign;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');

delta = 0.1; % 0.1 V bin size
bins = [-2 : delta : 11];

subset = indices_flattop_no_disrupt;

slice_data = v_loop(subset) .* ipsign(subset);

counts = hist(slice_data, bins);
bar_handle = bar(bins, counts/sum(counts)/delta, 'hist');
%set(bar_handle, 'EdgeColor', 'none');
%set(bar_handle, 'FaceColor', 'k');
set(gca, 'linewidth', 1.5, 'fontsize', 12)

xlim([-1, 10]);
xlabel('V_{loop} [V]', 'FontSize', 14);
ylabel('Probability [V^{-1}]', 'FontSize', 14);
title('C-Mod: Distribution of V_{loop}', 'FontSize', 16);

text(0.80, 0.89, 'All', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.82, 'non-disruptive', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.75, 'shots', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.68, '(flattop only)', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('cmod_vloop_no_disrupt_histogram.png', '-dpng');
