%db = set_database('d3drdb');
%retrieve_all_params

define_indices;

qq = indices_flattop_disrupt_in_flattop;

figure('menubar', 'none', 'toolbar', 'none', 'PaperPositionMode', 'auto');
plot(-smear(time_until_disrupt(qq), 0.0015), n_equal_1_mode(qq)*1.e3, '.');
xlim([-2.0, 0.0]); ylim([0,3]);
set(gca, 'fontsize', 12, 'linewidth', 1.5);
xlabel('Time until disrupt [s]', 'fontsize', 14);
ylabel('n=1 locked mode amplitude [mT]', 'fontsize', 14);
title('D3D locked mode amplitude vs time until disrupt', 'fontsize', 16);

print('d3d_locked_mode_vs_time_disrupt.png', '-dpng');
