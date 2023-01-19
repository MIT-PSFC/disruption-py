% <kmontes@mit.edu>
% MIT PSFC
% Created: Feb 2019

function [rho_diag,varargout] = efitmap_Rz_to_rho(diag_time,...
    diag_r,diag_z,pulse,varargin)
%EFITMAP_RZ_TO_RHO Convert absolute [R,z] channel coords to 1D midplane
%   normalized rho. Takes EFIT data from DIII-D MDS+, interpolates the 
%   user-specified positions onto flux surfaces, then maps to the midplane
%
%   Inputs:
%       - diag_time: timebase [s] for input diagnostic
%       - diag_r, diag_z: major radial and z coordinates [m] of measurement
%       - pulse: shot/discharge number on DIII-D
%       - type (optional): desired output. Is either normalized_rho
%       (default), absolute_rho (in meters), or major_radius (in meters)
%
%   Outputs:
%       - rho_diag: radial midplane coordinates after EFIT mapping at each
%       time slice in diag_time for each measuremente
%-------------------------------------------------------------------------

% Check input dimensions (time on left, channel number on right)
if size(diag_r,2)==1; diag_r = diag_r'; end
if size(diag_z,2)==1; diag_z = diag_z'; end
if size(diag_time,1)==1; diag_time = diag_time'; end

% Grab EFIT data for this pulse from MDS+
if nargin > 5
    trees = varargin{2};
    EFIT = load_efit(pulse,trees);
else
    EFIT = load_efit(pulse);
end

if ~isstruct(EFIT)
    disp('Failed to load EFIT data...')
    return
end

%% Find EFIT slices to use for each Thomson slice

slice_to_use = zeros(size(diag_time));
% If thomson starts before EFIT (often does), then use the first valid EFIT slice for early Thomson data.
% Bad, but nothing better. User should be able to notice this happening from a contour.
early_indices = find(diag_time < min(EFIT.time));
if length(early_indices) > 0
    slice_to_use(early_indices) = 1;
    first_ts = max(early_indices)+1;
else;
    first_ts = 1;
end
% If Thomson ends after EFIT (also often happens), then use the last valid EFIT slice for late Thomson data.
% Thomson after EFIT probably happened after the flattop and you don't care, anyway. It should be apparent
% that the same EFIT is being used if you look at a contour.
late_indices = find(diag_time > max(EFIT.time));
if length(late_indices) > 0
    slice_to_use(late_indices) = length(EFIT.time);
    last_ts = min(late_indices) - 1;
else;
    last_ts = length(diag_time);
end

% Acquire list of diag time slices w/in EFIT time range; should find 
% closest EFIT for each one
diag_slices = [first_ts:1:last_ts];
for i=diag_slices
    [~,closestIndex] = min(abs(EFIT.time-diag_time(i)));
    slice_to_use(i) = closestIndex;
end

%% Interpolation (This section translated from David Eldon's OMFITthomson map routine in Python)
% Okay, great. We should know which EFIT slice is best for each TS slice. You could interpolate the EFIT in
% time to be fancier, but this probably doesn't matter. It's fine. Let's go!
%
% We're going to do a fast spatial interpolation.
% This is optimized using specific information about the Thomson system: that almost all of the chords
% are in the same radial position. That means that we can calculate the coefficients for the R
% interpolation once for each unique R value (there is only 1 for the core, and only 1 for divertor.
% Tangential is different). We do the R interpolation once, get psi_N vs Z at the correct R, and then
% do a simple 1D interpolation. Fast! For the tangential system, it will loop through the 6 unique R
% values instead of doing just 1 unique R value for core & div. Oh well.

% Make an array to catch psin at the TS locations vs EFIT time base. This guy needs to be converted to the
% TS time base later. The mapping for this conversion was calculated earlier.
psin_diag_almost = zeros(length(EFIT.time),length(diag_z));

for r = unique(diag_r)
    dr = r - EFIT.r;
    right = find(EFIT.r > r,1); % Find closest EFIT R on the left & 
    left = right - 1;              % right of TS R value
    if EFIT.r(right) == r
        psin_slice = squeeze(EFIT.psin(:,right,:));
    else
        dlr = EFIT.r(right)-EFIT.r(left); % Find the spacing between the EFIT R values on L and R
        left_weight = -dr(right)/dlr; % Large diff between TS & right EFIT ==> large weight on left value
        right_weight = dr(left)/dlr; % Large diff between TS & left EFIT ==> large weight on right value
        % Weighted sum of psin from left and right sides:
        psin_slice = squeeze(EFIT.psin(:,left,:)*left_weight + ...
            EFIT.psin(:,right,:)*right_weight); % This reduces psin from psin(R,Z) to psin(Z) 
    end
    % Now narrow down the z-positions to iterate through for this r-value
    z_vals = diag_z(find(diag_r==r));
    for z = unique(z_vals)
        % Do the Z interpolation for each unique TS_z value; Same idea as
        % before, except it reduces a line to a point instead of a grid to a line.
        dz = z - EFIT.z;
        right = find(EFIT.z > z,1); left = right - 1;
        if EFIT.z(right) == z
            psin_pt = psin_slice(:,right);
        else
            dlr = EFIT.z(right) - EFIT.z(left);
            left_weight = -dz(right)/dlr;
            right_weight = dz(left)/dlr;
            psin_pt = psin_slice(:,left)*left_weight + ...
                psin_slice(:,right)*right_weight;
        end
        % Now we map this back into the shape of TS data. Remember, this interpolation operated on unique
        % values within the TS coordinates, not on each point. This result here does not necessarily have
        % compatible dimensions with the original TS, so we have to map it in.
        psin_diag_almost(:,find(diag_r==r & diag_z == z)) = psin_pt;
    end      
end

% Normalized flux on TS timebase for each channel
psin_diag = psin_diag_almost(slice_to_use,:);

% Finish the mapping by obtaining normalized rho at each TS time and channel
rhovn_diag_almost = interp1(EFIT.time,EFIT.rhovn,diag_time); % Get rhovn on TS timebase
rhovn_diag = NaN(size(psin_diag));
psin_interp = linspace(0,1,size(EFIT.rhovn,2)); % implied psin grid for rhovn 
for i = 1:length(diag_time) % Interpolate again to get rhovn on same psin base
    rhovn_diag(i,:) = interp1(psin_interp,rhovn_diag_almost(i,:),psin_diag(i,:));
end

%% Convert rhovn_diag to the desired output

% Write if/then cases here to handle different arguments 'type' (maybe
% normalized_rho, absolute_rho, etc...)

if nargin > 4
    type = varargin{1};
    if strcmp(type,'absolute_rho'); % not sure about this -> depends on how rhovn and aminor are defined
        aminor_diag = interp1(EFIT.aminor_time,EFIT.aminor,diag_time);
        rho_diag = aminor_diag.*rhovn_diag;
    elseif strcmp(type,'major_radius');
	[~,idx] = min(abs(EFIT.z)); % find magnetic axis
	midplane_flux = squeeze(EFIT.psin(:,:,idx)); % poloidal flux on R grid @ z=0 on EFIT timebase
	midplane_flux = interp1(EFIT.time,midplane_flux,diag_time); % interpolate onto TS timebase
	rmaxis_interp = interp1(EFIT.time,EFIT.rmaxis,diag_time); % interpolated R_0 onto TS timebase
	rho_diag = NaN(size(rhovn_diag));
	good_indices = find(~all(isnan(midplane_flux')));
	for i=1:length(diag_time);
	    if ismember(i,good_indices);
		low_field_side_indices = find(EFIT.r>rmaxis_interp(i));
		rho_diag(i,:) = interp1(midplane_flux(i,low_field_side_indices), ...
			EFIT.r(low_field_side_indices)',psin_diag(i,:)); % get R vals for diag
	    end	
	end
    else % outputs normalized rho_Vn, as stored in the EFIT tree
        rho_diag = rhovn_diag;
    end
else % outputs normalized rho_Vn, as stored in the EFIT tree
    rho_diag = rhovn_diag;
end

if nargout > 1
    varargout{1} = EFIT;
end
