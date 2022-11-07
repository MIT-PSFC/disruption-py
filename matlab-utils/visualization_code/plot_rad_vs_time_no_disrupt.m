%retrieve_all_params;
define_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');
plot(smear(time(indices_flattop_no_disrupt), 0.0),...
  p_rad(indices_flattop_no_disrupt)/1.e6, '.');
xlim([0, 2]); ylim([0,5]);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
xlabel('Time [s]', 'FontSize', 14);
ylabel('P_{rad} [MW]', 'FontSize', 14);
title('C-Mod: P_{rad}', 'FontSize', 16);

print('cmod_rad_vs_time_no_disrupt.png', '-dpng');
