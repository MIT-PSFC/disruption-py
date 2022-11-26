function [Te_peaking_factor, ne_peaking_factor, profile] = ...
    get_Thomson_peaking_factors(shot, timebase,pftype);

% This function retrieves data from the Thomson scattering diagnostic on 
% DIII-D and uses it to calculate peaking factors describing the shapes
% of the Te and ne profiles..
% 
% Dependent on the following routines:
%   - load_thomson.m
%   - efitmap_Rz_to_rho.m
%   - load_efit.m
%
% Inputs:
%   shot = shot number
%   timebase = times to calculate the powers
%   sigtype = either 'temp' or 'density' string (specifies which peaking
%   factor to calculate)
% Outputs:
%   peaking_factor = ratio of 'core' to 'edge' values (core_vs_edge) or
%   'core' to 'all' values (core_vs_all)
%   profile = structure containing raw and interpolated data for debugging
%
% Author: Kevin Montes (04-13-2019)

% Initialize all outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

Te_peaking_factor = NaN(size(timebase));
ne_peaking_factor = NaN(size(timebase));
profile = [];

%% Make sure inputs are given in desired format
if (size(timebase,1) > 1);
  timebase = timebase;
else;
  timebase = transpose(timebase);
end;

if ~strcmp(pftype,'core_vs_edge') && ~strcmp(pftype,'core_vs_all')
    disp([pftype,' is not a valid "pftype" input. ', ...
        'Need "core_vs_edge" or "core_vs_all"'])
    return
end

% EFIT runs that are done automatically after each shot are usually put in
% the EFIT01 and EFIT02 trees.  Additional EFIT runs using non-standard
% settings and/or time sampling are often done by individuals, and these go
% into trees named EFIT03, EFIT04, etcetera.  This includes specialized
% EFIT runs that we have done for use with our disruption warning database.
% An SQL database of available EFIT trees exists, from which our specific
% EFIT tree can be selected by means of the 'run_by' and/or 'runtype'
% keywords.  All of our EFIT runs have runtype = 'DIS' (for disruption).

efittrees  = select_efit_trees(shot, '', 'DIS');
%efittrees  = select_efit_trees(shot, 'granetzr', 'DIS');
%efittrees = select_efit_trees(shot, {'granetzr', 'reac'}, 'DIS');

if isempty(efittrees);
  fprintf(1, 'No disruption EFIT tree for this shot\n');
  return;
end;
tree=char(efittrees(end,:));

try
    TS = get_thomson(shot); % Load Thomson scattering data
catch ME
    disp(ME)
    return
end

% Map Thomson channel coordinates to the midplane
try
    [rho_TS,EFIT] = efitmap_Rz_to_rho(TS.time,...
        TS.R,TS.z,double(shot),'major_radius',{tree});
catch ME
    disp(ME)
    return
end

% Store TS measurements in a structure called 'data'
data.t = TS.time;
data.Te = TS.Te';
data.ne = TS.ne';
data.x = rho_TS';
data.xtime = TS.time;

% Get magnetic axis and minor radius from EFIT structure
RMAG.signal = EFIT.rmaxis;
RMAG.time = EFIT.time;
MINRAD.signal = EFIT.aminor;
MINRAD.time = EFIT.aminor_time;

%% Check data shape and map onto input timebase
if size(data.x,1) ~= size(data.Te,1)
    if size(data.x,1) == size(data.Te,2) && numel(data.Te) == size(data.Te,1)
        data.Te = data.Te';
	    data.ne = data.ne';
    else
        warning('Data incomplete!')
        return
    end
end

if ~isequal(data.t,timebase)
    data.Te = data.Te(:,dsearchn(data.t(:),timebase));
    data.ne = data.ne(:,dsearchn(data.t(:),timebase));
end

%% Output the raw data for debugging purposes
save_mode_data = 'off';
if strcmp(save_mode_data,'on')
    profile.data = data; % set 'on' only to download diagnostics full coordinates!
end

%% Radial coordinate interpolation 
% -------------------------------
if isfield(data,'xtime')
    data.xinterp = interp1(data.xtime,data.x',timebase,'linear','extrap')';
else
    data.xinterp = repmat(data.x,1,numel(timebase));
end

%% Get inner and outer radii
% Get magnetic axis (INNER) radius on the input timebase
if isfield(RMAG,'time')
    inn = interp1(RMAG.time,RMAG.signal,timebase,'linear','extrap')';
else
    inn = RMAG.signal;
end

% Get plasma minor (OUTER) radius on the input timebase
if isfield(MINRAD,'time')
    out = interp1(MINRAD.time,MINRAD.signal, timebase,'linear','extrap')';
else
    out = MINRAD.signal;
end

% Get major radial coordinate of the OUTER boundary (if not bolometer)
out = inn+out; % evaluated wrt the inner coordinate (magnetic axis radial position)

%% Account for margins (specified below, empty by default)
inner_marg = []; outer_marg = [];
% Get range of major radial coords and exclude margin (outer 30%?) from profile, if specified
prof_range = out - inn;
if ~isempty(inner_marg)
    inn = inn+inner_marg*prof_range; 
end
if ~isempty(outer_marg)
    out = out-outer_marg*prof_range; 
end
prof_range = out - inn;

%% General rule (find indices of radii lying within inner/outer bounds)
% ------------
idx_rad_coord = bsxfun(@ge,data.xinterp,inn) & bsxfun(@le,data.xinterp,out);

z = {'Te','ne'};
for i=1:2
    signal = [z{i},'_signal'];
    % Profile calculation
    % -------------------
    profile.(signal) = nan(size(data.(z{i}))); 
    profile.(signal)(idx_rad_coord) = data.(z{i})(idx_rad_coord); 
    profile.time = timebase;

    %% Peaking factor section 
    % ----------------------
    % Set core & edge regions (consistent w/ Ale's on other tokamaks)
    pc_core_def = 0.30; % constant for all profiles (~30% of minor radius)
    pc_edge_def = [0.8, 0.95]; % within 80% - 95% of minor radius

    % Core definition
    % ---------------
    idx_core = bsxfun(@ge, data.xinterp, inn) & ...
        bsxfun(@le, data.xinterp, inn+pc_core_def*prof_range);
    core_mat = nan(size(profile.(signal))); 
    core_mat(idx_core) = profile.(signal)(idx_core);

    % Edge definition
    % ---------------
    idx_edge = bsxfun(@ge, data.xinterp, inn+pc_edge_def(1)*prof_range) & ...
        bsxfun(@le, data.xinterp, inn+pc_edge_def(2)*prof_range);
    edge_mat = nan(size(profile.(signal)));  
    edge_mat(idx_edge) = profile.(signal)(idx_edge);  

    profile.(signal) = abs(profile.(signal));
    profile.idx_core = idx_core;
    profile.idx_edge = idx_edge;

    if strcmp(pftype,'core_vs_all')
        for j = 1:size(core_mat,2)
            peaking_factor.(z{i})(j) = mean(core_mat(~isnan(core_mat(:,j)),j))./mean(profile.(signal)(~isnan(profile.(signal)(:,j)),j));
            profile.([z{i},'_core'])(j) = mean(core_mat(~isnan(core_mat(:,j)),j));
        end
    elseif strcmp(pftype,'core_vs_edge')
        for j = 1:size(core_mat,2)
            peaking_factor.(z{i})(j) = mean(core_mat(~isnan(core_mat(:,j)),j))./mean(edge_mat(~isnan(edge_mat(:,j)),j));
            profile.([z{i},'_core'])(j) = mean(core_mat(~isnan(core_mat(:,j)),j));
        end
    end
end

% Assign outputs
Te_peaking_factor = peaking_factor.Te;
ne_peaking_factor = peaking_factor.ne;
