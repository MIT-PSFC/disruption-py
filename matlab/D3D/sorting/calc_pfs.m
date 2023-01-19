function [PF,TS,Prad,EFIT] = calc_pfs(pulse,data_source)
%CALC_PEAKING_FACTORS Retrieves DIII-D peaking factors for Te, ne, and Prad
%   Author: Kevin Montes
%   Date: March 2019
%	Inputs:
%		- pulse: shot number for which to pull data
%		- data_source: source from which to grab raw data ('mdsplus' or '.mat')
%----------------------------------------------------

%% Load raw data
%-------------------------
if strcmp(data_source,'mdsplus')
	% Pull all data from the DIII-D MDS+ server (requires DIII-D sign-in)
	disp('Loading Thomson data...')
	addpath(['/fusion/projects/disruption_warning/peaking_factors_d3d/' ...
				'Physics-based_indicators/DIAG_parameterization'])
	TS = load_ne_Te(pulse,'blessed');
	disp('Mapping Thomson data using EFIT...')
	efit_trees = {'EFIT01','EFIT02','EFIT03','EFIT04','EFIT05'}; % cell array of EFIT tree names to try for mapping
	[TS.rho,EFIT] = efitmap_Rz_to_rho(TS.time,...
		TS.R,TS.z,pulse,'normalized',efit_trees);
	disp('Loading bolometer data...')
	smoothing_window = 40; % window for filtering raw signal in [ms]
	Prad = load_Prad(pulse,'lower',smoothing_window,EFIT);
elseif strcmp(data_source,'.mat')
	% Pull all data from local .mat files
	load(['EFIT_',num2str(pulse)]) % structure holding needed measurements from EFIT
	load(['TS_',num2str(pulse)]) % structure holding Te, ne, and rho from Thomson
	load(['Prad_',num2str(pulse)]) % structure holding Prad, from bolometer
end

% Find core bin for Thomson and calculate Te, ne peaking factors
%-----------------------------------------------------------------
core_mask = TS.rho<0.3; % core channels are within 30% of normalized rho
Te_core = TS.Te; Te_core(~core_mask) = nan; 
ne_core = TS.ne; ne_core(~core_mask) = nan; % omit other chans from core by setting to NaN
PF.Te.pf = mean(Te_core,2,'omitnan')./mean(TS.Te,2,'omitnan');
PF.ne.pf = mean(ne_core,2,'omitnan')./mean(TS.ne,2,'omitnan'); % 'core vs all' peaking factors
PF.Te.time = TS.time;
PF.ne.time = TS.time; % record raw diagnostic time as the peaking factor time

% Find core bin for bolometer and calculate Prad CVA, X-DIV peaking factors
%--------------------------------------------------------------------------

% Interpolate zmaxis and channel intersects x onto the bolometer timebase
zmaxis = interp1(EFIT.time,EFIT.zmaxis,Prad.t);
zmaxis = repmat(zmaxis',1,size(Prad.x,1))';
Prad.xinterp = interp1(Prad.xtime,Prad.x',Prad.t,'linear','extrap')';

% Determine the bolometer channels falling in the 'core' bin
vert_range = 3; % vertical range of the DIII-D cross section in meters
Prad_core_def = 0.06; % percentage of cross-section defining the core margin
core_bin = (Prad.xinterp < zmaxis + Prad_core_def*vert_range) & ...
				(Prad.xinterp > zmaxis - Prad_core_def*vert_range);

% Designate the divertor bin and find all 'other' channels not in that bin
div_channels = [1:7]+24; % lower bolometer fan channel numbers covering divertor
div_indices = find(ismember(Prad.ch_avail,div_channels));
other_indices = find(~ismember([1:size(Prad.x,1)],div_indices));

% Grab Prad measurements for each needed set of channels
Prad_core = Prad.z; Prad_core(~core_bin) = NaN; % core (P_j for j in C)
Prad_all_but_core = Prad.z; Prad_all_but_core(core_bin) = NaN; % All P_j excluding j in C
Prad_div = Prad.z(div_indices,:); % divertor (P_j for j in D)
Prad_all_but_div = Prad.z(other_indices,:); % All P_j excluding j in D

% Calculate Prad peaking factors
PF.Prad_CVA.pf = (mean(Prad_core,1,'omitnan')./mean(Prad_all_but_div,1,'omitnan'))';
PF.Prad_XDIV.pf = (mean(Prad_div,1,'omitnan')./mean(Prad_all_but_core,1,'omitnan'))';
PF.Prad_CVA.time = Prad.t';
PF.Prad_XDIV.time = Prad.t';
