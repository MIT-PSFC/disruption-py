%retrieve_all_params;
define_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');
plot(smear(time(indices_flattop_no_disrupt), 0.0),...
  radiated_fraction(indices_flattop_no_disrupt), '.');
xlim([0, 2]); ylim([0,2]);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
xlabel('Time [s]', 'FontSize', 14);
ylabel('P_{rad}/P_{input}', 'FontSize', 14);
title('C-Mod: P_{rad}/P_{input}', 'FontSize', 16);

print('cmod_rad_frac_vs_time_no_disrupt.png', '-dpng');
