% <kmontes@mit.edu>
% MIT PSFC
% Created: Feb 2020

function [data,varargout] = efit_Rz_interp(data,EFIT,varargin)
%EFITMAP_RZ_TO_RHO Convert absolute [R,z] coordinates to either normalized flux
% values, rhoVn, or major radius values. Takes EFIT data from EAST MDS+ and
% interpolates the user-specified positions onto flux surfaces
%
%   Inputs:
%       - data: a structure holding data of interest with the following fields
%				- data.R: a 1D vector or 2D matrix (time dependent) of radial coordinates
%				- data.z: a 1D vector or scalar value (if constant) of vertical position coordinates
%				- data.time: a 1D vector of times at which data is sampled
%       - EFIT: a structure returned by 'load_efit_east.m' routine, corresponding
%				to the shot for which the 'data' structure is loaded 
%       - type (optional): desired output. Is either normalized flux 'psin'
%       (default), 'rhoVn', or 'major_radius' (in meters)
%
%   Outputs:
%       - data: returns the original data structure, now with an additional
%				field called 'psin', 'rhoVn', or 'major_radius' (depending on
%				the specified interpolation type)
%-------------------------------------------------------------------------

% Check input dimensions (time on left, channel number on right)
if isvector(data.R)
	if size(data.R,2)==1; data.R = data.R'; end
	if size(data.z,2)==1; data.z = data.z'; end
	T = repmat(data.time,1,length(data.R)); % time of full grid
	R = repmat(data.R,length(data.time),1); % major radius of full grid
	Z = repmat(data.z,length(data.time),1); % vertical position of full grid
elseif ismatrix(data.R)
	if size(data.R,1)~=length(data.time) & size(data.R,2)==length(data.time)
		data.R = data.R'; 
	end
	data.z = unique(data.z);
	if length(data.z)~=1
		disp('Error: Expected 1 z value for reflectometer ...')
		return
	end
	T = repmat(data.time,1,size(data.R,2)); % time of full grid
	R = data.R;
	Z = repmat(data.z,length(data.time),size(data.R,2)); % vertical position of full grid
end

if size(data.time,1)==1; data.time = data.time'; end

% Implement a 3D (time,radial,vertical) gridded interpolation
gridVecs = {EFIT.time,EFIT.r,EFIT.z};
F = griddedInterpolant(gridVecs,EFIT.psin,'linear','none'); % interpolant function

% Apply interpolant to diagnostic data and return outputs as a new structure field
psin = F(T,R,Z);
data.psin = psin;

% Get rhovn using the interpolant stored in EFIT, then save this as another field in 'data'
rhovn_diag_almost = interp1(EFIT.time,EFIT.rhovn,data.time); % Get rhovn on TS timebase
rhovn_diag = NaN(size(psin));
psin_interp = linspace(0,1,size(EFIT.rhovn,2)); % implied psin grid for rhovn 
for i = 1:size(psin,1) % Interpolate again to get rhovn on same psin base
    rhovn_diag(i,:) = interp1(psin_interp,rhovn_diag_almost(i,:),psin(i,:));
end
data.rhovn = rhovn_diag;

if nargout > 1
    varargout{1} = EFIT;
end
