function [Te_PF, ne_PF, Rad_CVA, Rad_XDIV] = get_peaking_factors_d3d(shot,timebase)
%_______________________________________________________________________
% This function calculates peaking factors for the shot number
% given by the user corresponding to the times in the given timebase.
% Electron temperature (Te_PF) and density (ne_PF) profile peaking 
% factors are taken from Thomson scattering measurements, and the peaking
% factors describing radiated power distributions (Rad_CVA and Rad_XDIV)
% are taken from the 2pi foil bolometer system. 
%
% Inputs:
%	- shot: shot number (integer or double type)
%	- timebase: a 1d vector of (sorted) times at which to calculate PFs 
% Outputs:
%	- Te_PF: core Te vs. all Te peaking factor from Thomson scattering
%	- ne_PF: core ne vs. all ne peaking factor from Thomson scattering
%	- Rad_CVA: core vs. all brightness/P_rad metric from bolometer
%	- Rad_XDIV: divertor vs. all brightness/P_rad metric from bolometer
%
% Author: Kevin Montes
% Date: May 19, 2020
%_______________________________________________________________________

% Specify peaking factor variables
%-------------------------------------------
TS_data_type = 'blessed'; % either 'blessed', 'unblessed', or 'ptdata'
TS_radius = 'rhovn'; % metric to use for core/edge binning (either 'psin' or 'rhovn')
TS_core_margin = 0.3; % TS_radius value defining boundary of 'core' region (between 0 and 1)
TS_radial_range = [0,1]; % All data outside this range excluded. For example, psin=0 at magnetic axis and 1 at separatrix.
TS_equispaced = false; % set to true to interpolate TS_channel data onto equispaced radial grid
bolometer_fan = 'custom'; % fan to use for P_rad peaking factors (either 'lower', 'upper', or 'custom')
div_channels = [3:7]+24; % array of bolometer fan channel numbers covering divertor (upper fan: 1->24, lower fan: 25:48) 
smoothing_window = 40; % time window for filtering raw bolometer signal in [ms]
Prad_core_def = 0.06; % percentage of DIII-D veritcal extent defining the core margin
Prad_metric = 'brightness';%'brightness'; % either 'brightness' or 'power' ('z')

% Add libraries & paths to necessary scripts
% ------------------------------------------
addpath('/fusion/projects/disruption_warning/software/peaking_factors_d3d/Physics-based_indicators/DIAG_parameterization'); % path [1]
addpath('/fusion/projects/disruption_warning/software/peaking_factors_d3d/shared_scripts'); % path [2]

% Make sure input timebase is column vector, and input shot number is double type 
%-------------------------------------------------------------------------------
if isvector(timebase)
	Te_PF = NaN(size(timebase));	%--------------------------------------------------
	ne_PF = NaN(size(timebase));	% Initialize each peaking factor to NaN values. This
	Rad_CVA = NaN(size(timebase));  % is what will be returned if the routine fails.
	Rad_XDIV = NaN(size(timebase)); %--------------------------------------------------
	original_format = 'column';
	if ~iscolumn(timebase); 
		timebase = timebase';
		original_format = 'row';
	end
else
	disp("ERROR: Input 'timebase' must be a 1D vector"); return
end
if ~isa(shot,'double')
    shot = cast(shot,'double');
end

% Select disruption EFIT tree and retrieve a structure with relevant variables
%-----------------------------------------------------------------------------
% EFIT runs that are done automatically after each shot are usually put in
% the EFIT01 and EFIT02 trees.  Additional EFIT runs using non-standard
% settings and/or time sampling are often done by individuals, and these go
% into trees named EFIT03, EFIT04, etcetera.  This includes specialized
% EFIT runs that we have done for use with our disruption warning database.
% An SQL database of available EFIT trees exists, from which our specific
% EFIT tree can be selected by means of the 'run_by' and/or 'runtype'
% keywords.  All of our EFIT runs have runtype = 'DIS' (for disruption).

efittrees  = select_efit_trees(shot, '', 'DIS');
if isempty(efittrees);
  fprintf(1, 'No disruption EFIT tree for this shot\n');
  return;
end;
EFIT = load_efit(shot,efittrees); % see path [1]
if ~isstruct(EFIT)
	disp(['ERROR: Failed to load EFIT structure for shot #' num2str(shot)])
	return
end

% Get raw Te, ne, and P_rad data (see paths [1] and [2])
%------------------------------------------------------------------------
try
	TS = load_ne_Te(shot,TS_data_type); % load Thomson scattering data
	TS = efit_Rz_interp(TS,EFIT); % map Thomson (R,Z) measurements to flux (psi) values
catch
	TS = 0;
end
try
	Prad = load_Prad(shot,bolometer_fan,smoothing_window,EFIT); % load bolometer data
catch
	Prad = 0;
end
if ~isstruct(TS) & ~isstruct(Prad)
	disp(['ERROR: Both TS and bolometer data missing for shot #' num2str(shot)])
	return
end

%________________________________
% Te & ne Peaking Factors
%________________________________
if isstruct(TS)
	try
		% Drop data outside of valid range
		%----------------------------------------------
		out_of_bounds = (TS.(TS_radius)<TS_radial_range(1) & TS.(TS_radius)>TS_radial_range(2));
		TS.Te(out_of_bounds) = nan; TS.ne(out_of_bounds) = nan;
		TS.Te(isnan(TS.(TS_radius))) = nan;
		TS.ne(isnan(TS.(TS_radius))) = nan;

		% Interpolate onto uniform radial base if needed
		%-----------------------------------------------------------------
		if TS_equispaced
			for i=1:length(TS.time)
				no_nans = find(~isnan(TS.Te(i,:)) & ~isnan(TS.ne(i,:)));	
				if ~isempty(no_nans) & length(no_nans)>1
					radii = TS.(TS_radius)(i,no_nans);
					if length(radii)>2
						rad_coord_interp = linspace(min(radii),max(radii),numel(radii));
						TS.Te(i,no_nans) = interp1(radii,TS.Te(i,no_nans),rad_coord_interp,'pchip');
						TS.ne(i,no_nans) = interp1(radii,TS.ne(i,no_nans),rad_coord_interp,'pchip');
						TS.(TS_radius)(i,no_nans) = rad_coord_interp;
					end
				end
			end
		end

		% Find core bin for Thomson and calculate Te, ne peaking factors
		%-----------------------------------------------------------------
		core_mask = TS.(TS_radius)<TS_core_margin ; % mask array that is true only in the core region
		Te_core = TS.Te; Te_core(~core_mask) = nan; 
		ne_core = TS.ne; ne_core(~core_mask) = nan; % omit other chans from core by setting to NaN
		Te_PF = mean(Te_core,2,'omitnan')./mean(TS.Te,2,'omitnan');
		ne_PF = mean(ne_core,2,'omitnan')./mean(TS.ne,2,'omitnan'); % 'core vs all' peaking factors
		Te_PF = interp1(TS.time,Te_PF,timebase);
		ne_PF = interp1(TS.time,ne_PF,timebase); % interpolate onto given timebase

		% Switch format back if needed
		if strcmp(original_format,'row'); Te_PF = Te_PF'; ne_PF = ne_PF'; end
	catch
		disp(['ERROR: Failed to calculate Te, ne peaking factors for shot #' num2str(shot)])
		return
	end
end

%________________________________
% Prad CVA, X-DIV Peaking Factors
%________________________________
if isstruct(Prad)
	try
		% Interpolate zmaxis and channel intersects x onto the bolometer timebase
		zmaxis = interp1(EFIT.time,EFIT.zmaxis,Prad.t);
		zmaxis = repmat(zmaxis',1,size(Prad.x,1))';
		Prad.xinterp = interp1(Prad.xtime,Prad.x',Prad.t,'linear','extrap')';

		% Determine the bolometer channels falling in the 'core' bin
		vert_range = 3; % vertical range of the DIII-D cross section in meters
		core_bin = (Prad.xinterp < zmaxis + Prad_core_def*vert_range) & ...
						(Prad.xinterp > zmaxis - Prad_core_def*vert_range);

		% Designate the divertor bin and find all 'other' channels not in that bin
		div_indices = find(ismember(Prad.ch_avail,div_channels));
		other_indices = find(~ismember([1:size(Prad.x,1)],div_indices));

		% Grab Prad measurements for each needed set of channels
		Prad_core = Prad.(Prad_metric); Prad_core(~core_bin) = NaN; % core (P_j for j in C)
		Prad_all_but_core = Prad.(Prad_metric); Prad_all_but_core(core_bin) = NaN; % All P_j excluding j in C
		Prad_div = Prad.(Prad_metric)(div_indices,:); % divertor (P_j for j in D)
		Prad_all_but_div = Prad.(Prad_metric)(other_indices,:); % All P_j excluding j in D

		% Calculate Prad peaking factors
		Rad_CVA = (mean(Prad_core,1,'omitnan')./mean(Prad_all_but_div,1,'omitnan'))';
		Rad_XDIV = (mean(Prad_div,1,'omitnan')./mean(Prad_all_but_core,1,'omitnan'))';
		Rad_CVA = interp1(Prad.t',Rad_CVA,timebase);
		Rad_XDIV = interp1(Prad.t',Rad_XDIV,timebase);

		% Switch format back if needed
		if strcmp(original_format,'row'); Rad_CVA = Rad_CVA'; Rad_XDIV = Rad_XDIV'; end
	catch
		disp(['ERROR: Failed to calculate P_rad peaking factors for shot #' num2str(shot)])
		return
	end
end
%Rad_CVA = Prad;
