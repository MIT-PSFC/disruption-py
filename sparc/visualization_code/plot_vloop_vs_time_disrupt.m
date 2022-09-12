%retrieve_all_params;
define_indices;
get_ip_sign;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');
plot([-.025,0], [0,0], 'k', 'linewidth', 1.5);
xlim([-.025,0]); ylim([-1,10]);
hold on;
plot(-smear(time_until_disrupt(indices_flattop_disrupt_in_flattop),.0005),...
  v_loop(indices_flattop_disrupt_in_flattop) .* ...
  ipsign(indices_flattop_disrupt_in_flattop), '.');
set(gca, 'linewidth', 1.5, 'fontsize', 12)
xlabel('Time before disrupt [s]', 'FontSize', 14);
ylabel('Loop voltage [V]', 'FontSize', 14);
title('C-Mod: Loop voltage', 'FontSize', 16);

print('cmod_vloop_vs_time_disrupt.png', '-dpng');
