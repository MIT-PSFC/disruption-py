%{
retrieve_all_params;
define_indices;
%}

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-smear(time_until_disrupt(indices_flattop_disrupt_in_flattop))*1.e3, ...
  li(indices_flattop_disrupt_in_flattop),'.');
xlim([-50,0]);
ylim([0,2]);
set(gca, 'fontsize', 11);
xlabel('Time before disrupt [ms]', 'FontSize', 14);
ylabel('li (from EFIT)', 'FontSize', 14);
title('li vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng plot_li.png;
