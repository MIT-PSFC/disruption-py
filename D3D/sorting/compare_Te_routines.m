shot = 184250;
efittree  = select_efit_trees(shot, '', 'DIS');
t_plot = 2; % [s]
timebase = [0.67:0.005:5.12];
[~,t_indx] = min(abs(timebase-t_plot));

% Get legacy profile information
Diag_array = {'TS_{TE}';'TS_{NE}'};
Profiles = cell(1,numel(Diag_array));
Profiles(1,:) = DIS_tool_profiles_gen(shot,timebase,Diag_array,'lower',efittree);
data = Profiles{1}{2}.data;
signal = Profiles{1}{2}.signal;
idx_core = Profiles{1}{2}.idx_core;
core_mat = nan(size(signal));
core_mat(idx_core) = signal(idx_core);
core_val = mean(core_mat(~isnan(core_mat(:,t_indx)),t_indx));
all_val = mean(signal(~isnan(signal(:,t_indx)),t_indx))

% Get data using both routines
TS1 = load_ne_Te(shot,'blessed');
EFIT1 = load_efit(shot,efittree);
TS1 = efit_Rz_interp(TS1,EFIT);
TS2 = load_thomson(shot,'temp');
[rho2,EFIT2] = efitmap_Rz_to_rho(TS2.time,TS2.R,TS2.z,shot,'normalized',efittree);

% Get new profile information
[~,t_indx_new] = min(abs(TS1.time-t_plot));
core_mask = TS1.rhovn<0.3; % mask array that is true only in the core region
Te_core = TS1.Te; Te_core(~core_mask) = nan;
TS1.Te(isnan(TS1.rhovn)) = nan; 
new_core_val = mean(Te_core(t_indx_new,:),2,'omitnan');
new_all_val = mean(TS1.Te(t_indx_new,:),2,'omitnan')

% Plot both profiles for a specific time
[~,idx1] = min(abs(TS1.time-t_plot));
[~,idx2] = min(abs(TS2.time-t_plot));
figure(); hold on;
plot(TS1.rhovn(idx1,:),TS1.Te(idx1,:),'b.-','DisplayName','load\_ne\_Te')
plot(rho2(idx2,:),TS2.data(idx2,:),'rs','DisplayName','load\_thomson')
plot(xlim,[1,1]*core_val,'r-','DisplayName','legacy core')
plot(xlim,[1,1]*all_val,'r--','DisplayName','legacy all')
plot(xlim,[1,1]*new_core_val,'b-','DisplayName','new core')
plot(xlim,[1,1]*new_all_val,'b--','DisplayName','new all')
xlabel('rho')
ylabel('Te')
title(['Shot #' num2str(shot) ', t = ' num2str(t_plot) ' s'])
lgd = legend();
