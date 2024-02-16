%db = set_database('d3drdb');
%retrieve_all_params

define_indices;

qq = indices_flattop_no_disrupt;

figure('menubar', 'none', 'toolbar', 'none', 'PaperPositionMode', 'auto');
plot(time(qq), n_equal_1_mode(qq)*1.e3, '.');
xlim([0.0, 7.0]); ylim([0,3]);
set(gca, 'fontsize', 12, 'linewidth', 1.5);
xlabel('Time [s]', 'fontsize', 14);
ylabel('n=1 locked mode amplitude [mT]', 'fontsize', 14);
title('D3D locked mode amplitude (non-disruptions only)', 'fontsize', 16);

print('d3d_locked_mode_vs_time_no_disrupt.png', '-dpng');
