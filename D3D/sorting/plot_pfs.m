function [PF,TS,Prad,EFIT] = plot_pfs(pulse,data_source)
%CALC_PEAKING_FACTORS Retrieves DIII-D peaking factors for Te, ne, and Prad
%   Author: Kevin Montes
%   Date: March 2019
%	Inputs:
%		- pulse: shot number for which to pull data
%		- data_source: source from which to grab raw data ('realtime' or 'mdsplus')
%----------------------------------------------------

% Grab data
[PF,TS,Prad,EFIT] = calc_pfs(pulse,data_source);
min_time = 1.8;%min(EFIT.time);
max_time = 3.1;%max(EFIT.time);

% Generate figure
figure(); 
axs = []; hold off;

% Plot the peaking factors
ax = subplot(4,1,1);
hold on;
PF_names = fieldnames(PF);
PF_labels = {'$T_e$','$n_e$','$P_{rad}$ CVA','$P_{rad}$ XDIV'};
for i=1:length(PF_names)
	indices = intersect(find(PF.(PF_names{i}).time<=max_time),find(PF.(PF_names{i}).time>=min_time));
	plot(PF.(PF_names{i}).time(indices),PF.(PF_names{i}).pf(indices),'.-','DisplayName',PF_labels{i})
end
lgd = legend(); lgd.('Interpreter') = 'latex'; 
lgd.('FontSize')=14; %lgd.('NumColumns') = 2; 
hold off;
ylabel('Peaking Factors','FontSize',18)
set(ax,'XTickLabel',{''});
xlim([min_time,max_time])
axs = [axs,ax];

% Plot the density profiles on the 2nd panel
ii = find(TS.time<=max_time & TS.time>=min_time);
time_grid = repmat(TS.time',size(TS.rho,2),1);
ax = subplot(4,1,2);
ax.FontSize = 18;
surf(time_grid(:,ii),TS.rho(ii,:)',TS.ne(ii,:)'/1e19,'EdgeColor','None','facecolor','interp');
view([0 90])
ylabel('$\rho$','Interpreter','latex','FontSize',18)
set(ax,'XTickLabel',{''})
ylim([0,1])
axs = [axs,ax];
linkaxes(axs,'x')

% Plot the temperature profiles on the 3rd panel
ax = subplot(4,1,3);
ax.FontSize = 18;
surf(time_grid(:,ii),TS.rho(ii,:)',TS.Te(ii,:)'/1e3,'EdgeColor','None','facecolor','interp');
view([0 90])
ylabel('$\rho$','Interpreter','latex','FontSize',18)
set(ax,'XTickLabel',{''})
ylim([0,1])
axs = [axs,ax];
linkaxes(axs,'x')

% Plot the Prad distribution
ax = subplot(4,1,4);
ax.FontSize = 18;
indices = find(Prad.t>min_time & Prad.t<max_time);
imagesc(Prad.t(indices),[1:24],flipud(Prad.brightness(:,indices))/1e6,[0 1])
yticks([5 10 15 20])
yticklabels({'20','15','10','5'})
xlabel('Time [s]','FontSize',18)
ylabel('Channel #','FontSize',18)
axs = [axs,ax];
linkaxes(axs,'x')

% Format positions
% Adjust margins
margin_scale = 1.5;
marg = (0.5482-.3291-.1577)/2;
for i=1:length(axs)
	ax = axs(i); P = ax.Position; 
	ax.Position = [ax.Position(1)-.05,ax.Position(2)-marg,ax.Position(3),ax.Position(4)+marg];
end

%Changing cycle time for CPU16 from 500 microseconds to 0

% Set colormaps
%colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet')
x = axs(2).Position(1)+axs(2).Position(3)+0.01;
Pos = [x axs(2).Position(2) 0.03, axs(2).Position(4)];
c2 = colorbar(axs(2),'location','manual','Position',Pos);
c2.Label.String = '$n_e$ [$10^{19}$ $m^{-3}$]';
c2.Label.FontSize = 16;
c2.Label.Interpreter = 'latex'
Pos = [x axs(3).Position(2) 0.03, axs(3).Position(4)];
c3 = colorbar(axs(3),'location','manual','Position',Pos);
c3.Label.String = '$T_e$ [keV]';
c3.Label.FontSize = 16;
c3.Label.Interpreter = 'latex'
Pos = [x axs(4).Position(2) 0.03, axs(4).Position(4)];
c4 = colorbar(axs(4),'location','manual','Position',Pos);
c4.Label.String = 'Brightness [MW/$m^{-2}$]';
c4.Label.FontSize = 16;
c4.Label.Interpreter = 'latex'
%suptitle(['Shot #',num2str(pulse)])
