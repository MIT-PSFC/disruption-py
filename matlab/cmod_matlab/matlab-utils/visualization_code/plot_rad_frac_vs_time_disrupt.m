%retrieve_all_params;
define_indices;
get_ip_sign;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');
plot(-smear(time_until_disrupt(indices_flattop_disrupt_in_flattop),.0005),...
  radiated_fraction(indices_flattop_disrupt_in_flattop), '.');
xlim([-.025,0]); ylim([0,2]);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
xlabel('Time before disrupt [s]', 'FontSize', 14);
ylabel('P_{rad}/P_{input}', 'FontSize', 14);
title('C-Mod: P_{rad}/P_{input}', 'FontSize', 16);

print('cmod_rad_frac_vs_time_disrupt.png', '-dpng');
