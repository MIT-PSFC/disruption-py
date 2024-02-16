%retrieve_all_params;
define_indices;
get_ip_sign;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');
plot([0,2], [0,0], 'k', 'linewidth', 1.5);
xlim([0,2]); ylim([-1,10]);
hold on;
plot(smear(time(indices_flattop_no_disrupt),0.0),...
  v_loop(indices_flattop_no_disrupt) .* ...
  ipsign(indices_flattop_no_disrupt), '.');
set(gca, 'linewidth', 1.5, 'fontsize', 12)
xlabel('Time [s]', 'FontSize', 14);
ylabel('Loop voltage [V]', 'FontSize', 14);
title('C-Mod: Loop voltage', 'FontSize', 16);

print('cmod_vloop_vs_time_no_disrupt.png', '-dpng');
