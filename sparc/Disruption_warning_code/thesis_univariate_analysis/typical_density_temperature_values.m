%% Query database and get parameters

addpath('/home/montes/EPFL_collaboration/cmod_analysis_scripts');
db = set_database('logbook');
result = fetch(db,['select shot, time, time_until_disrupt, dipprog_dt, n_e ' ...
	'from disruption_warning order by shot, time'],'DataReturnFormat','cellarray');
shot = cell2mat(result(:,1));
time = cell2mat(result(:,2));
time_until_disrupt = cell2mat(result(:,3));
dipprog_dt = cell2mat(result(:,4));
n_e = cell2mat(result(:,5));
define_indices;

%% Plot density distribution
figure()
indices = intersect(indices_flattop,find(n_e>0 & n_e<1e21));
histogram(n_e(indices)/1e19,'Normalization','probability')
xlabel('$n_e\ [10^{19}]$','Interpreter','latex');
ylabel('probability')

%% Get flattop midpoint for each shot
shotlist = unique(shot(indices_flattop));
N = 1000;
ii = randsample(length(shotlist),N);
shot_subset = shotlist(ii);
Te0 = NaN(size(shot_subset));
mdsconnect('alcdata.psfc.mit.edu'); % Connect to C-Mod MDSplus server
for i=1:N
	shotnum = shot_subset(i);
	jj = intersect(indices_flattop,find(shot==shotnum));
	endpoints = [min(time(jj));max(time(jj))];
	t_mid = endpoints(1) + diff(endpoints)/2;
	disp(compose('Shot #%d : t_mid = %.2f s ... (%d/%d)',[shotnum,t_mid,i,N]))
	%disp(['Shot #' num2str(shotnum) ': t_mid = ' num2str(t_mid) ' s ...']);
	try
		TS = load_Thomson_cmod(shotnum);
		[~,closest] = min(abs(TS.time-t_mid));
		core_indx = find(TS.rho(closest,:)<0.3);
		Te0(i) = mean(TS.Te(closest,core_indx));
	catch ME
		disp(ME)
	end
end

%% Plot temperatures
disp(compose('Shots succeeded: %d',length(find(~isnan(Te0)))))
figure();
histogram(Te0);
xlabel('Te0')
