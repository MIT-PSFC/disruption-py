%% Pick a shot and time to plot profiles
shot = 161238;
t_plot = 3.5; % [s]
limits = [1,3];
addpath('/fusion/projects/disruption_warning/matlab_programs')
efittree  = select_efit_trees(shot, '', 'DIS');
EFIT = load_efit(shot,efittree);

% Get & plot blessed profile information
figure(); hold on;
source = {'ptdata','blessed','unblessed'};
for i=1:length(source)
	TS.(source{i}) = load_ne_Te(shot,source{i}); % retrieve data
	TS.(source{i}) = efit_Rz_interp(TS.(source{i}),EFIT); % map to rhovn
	[~,t_indx] = min(abs(TS.(source{i}).time-t_plot));
	plot(TS.(source{i}).rhovn(t_indx,:),TS.(source{i}).Te(t_indx,:),'.-','DisplayName',source{i})
	ii = find(TS.(source{i}).time>limits(1) & TS.(source{i}).time<limits(2));
	disp(['Number of time samples within limits, ' source{i} ': ' num2str(length(ii))])
end

xlabel('rho')
ylabel('Te')
title(['Shot #' num2str(shot) ', t = ' num2str(t_plot) ' s'])
lgd = legend();


% Plot all realtime profiles on separate figure
ii = find(TS.ptdata.time>limits(1) & TS.ptdata.time<limits(2));
figure();
mesh(TS.ptdata.rhovn,repmat(TS.ptdata.time,1,length(TS.ptdata.R)),TS.ptdata.Te)
xlabel('rho')
ylabel('time')
zlabel('Te')
%disp(['Number of times: ' num2str(length(ii))])

%{
% Get new profile information
[~,t_indx_new] = min(abs(TS1.time-t_plot));
core_mask = TS1.rhovn<0.3; % mask array that is true only in the core region
Te_core = TS1.Te; Te_core(~core_mask) = nan;
TS1.Te(isnan(TS1.rhovn)) = nan; 
new_core_val = mean(Te_core(t_indx_new,:),2,'omitnan');
new_all_val = mean(TS1.Te(t_indx_new,:),2,'omitnan')

% Plot both profiles for a specific time
plot(xlim,[1,1]*core_val,'r-','DisplayName','legacy core')
plot(xlim,[1,1]*all_val,'r--','DisplayName','legacy all')
plot(xlim,[1,1]*new_core_val,'b-','DisplayName','new core')
plot(xlim,[1,1]*new_all_val,'b--','DisplayName','new all')
%}
