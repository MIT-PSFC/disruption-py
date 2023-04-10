if (exist('time','var')==0 || exist('time_until_disrupt','var')==0 || ...
    exist('p_rad','var')==0 || exist('radiated_fraction','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

define_indices;

qq = indices_flattop_no_disrupt;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
set(gcf, 'PaperPositionMode', 'auto');
plot(smear(time(qq), 0.025), p_rad(qq)/1.e6,'.');
xlim([0,10]);
ylim([0,15]);
set(gca, 'fontsize', 12)
xlabel('Time [s]', 'FontSize', 14);
ylabel('P_{rad} [MW]', 'FontSize', 14);
title('P_{rad} vs time for non-disruptions (flattop ONLY)', 'FontSize', 16);
 
print -dpng d3d_prad_vs_time_no_disrupt.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
set(gcf, 'PaperPositionMode', 'auto');
plot(smear(time(qq), 0.025), radiated_fraction(qq),'.');
xlim([0,10]);
ylim([0,2]);
set(gca, 'fontsize', 12)
xlabel('Time [s]', 'FontSize', 14);
ylabel('P_{rad}/P_{input}', 'FontSize', 14);
title('P_{rad}/P_{input} vs time for non-disruptions (flattop ONLY)', ...
  'FontSize', 16);
 
print -dpng d3d_prad_frac_vs_time_no_disrupt.png;
