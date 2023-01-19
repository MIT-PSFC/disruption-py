if (exist('radiated_fraction','var')==0 || ...
    exist('time_until_disrupt','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

radiated_fraction_slow = p_rad_slow ./ (p_oh + p_lh + p_icrf);

get_flattop_indices;

no_disrupt_flattop = ...
  intersect(indices_flattop, find(isnan(time_until_disrupt)));

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-smear(time_until_disrupt(indices_flattop))*1.e3, ...         
  radiated_fraction_slow(indices_flattop),'.');
xlim([-50,0]);
ylim([0,2]);
set(gca, 'fontsize', 12);
xlabel('Time before disrupt [ms]', 'FontSize', 14);
ylabel('Radiated power fraction', 'FontSize', 14);
title('Radiated power fraction vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng pradfrac_slow_plot_3.png;
