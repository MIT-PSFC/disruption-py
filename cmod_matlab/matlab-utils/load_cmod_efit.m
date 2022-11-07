function EFIT = load_cmod_efit(pulse,varargin)
%LOAD_EFIT This function loads relevant EFIT data from a particular shot
%   Searches through the EFIT trees listed in efit_trees and outputs all
%   data in a MATLAB structure
% INPUTS:
%	- pulse: shot number
%	- efit_trees (optional): cell array of names of efit trees to 
%				 try, in order of priority
% OUTPUTS:
% 	- EFIT: a structure containing several quanitities pulled from
%		the efit files for the specified shot; if files fail
%		to load from MDSplus, EFIT = 0 is returned
%
% Author: Kevin Montes (04/2019)
% Email: kmontes@mit.edu
%----------------------------------------------------------------------

mdsconnect('alcdata.psfc.mit.edu'); % Start MDS+ connection for DIII-D

% List EFIT trees to try, in order of priority
if nargin > 1
    efit_trees = varargin{1};
else
    efit_trees = {'efit18','analysis'};
end

%% Pull the relevant data from MDS+
g_path = '\efit_geqdsk:'; %'\top.efit.results.geqdsk:';
a_path = '\efit_aeqdsk:'; %'\top.efit.results.aeqdsk:';
i = 0; loaded = false;
while i < length(efit_trees)
    i = i+1;
    % Open EFIT tree for this shot and grab the grid
    [shotopened, status]=mdsopen(efit_trees{i}, pulse);
    if mod(status,2)==1;
        % First grab poloidal flux time. All other measurements depend on
        % psi, so if this isn't retrievable we give up on this tree
        [EFIT.time,status] = mdsvalue(['dim_of(',g_path,'psirz,2)']);
        if mod(status,2)==1;
            EFIT.time = EFIT.time; % [s]
            EFIT.z = mdsvalue([g_path,'z']); % [m]
            EFIT.r = mdsvalue([g_path,'r']); % [m]
            EFIT.rhovn = mdsvalue([g_path,'rhovn']); % normalized r def 1
            psirz = mdsvalue([g_path,'psirz']); % measured poloidal flux
            ssimag = mdsvalue([g_path,'ssimag']);
            ssibry = mdsvalue([g_path,'ssibry']);
            EFIT.aminor = mdsvalue([a_path,'aminor']); % minor radius [m]
            EFIT.a_time = mdsvalue([a_path,'time']); % [s]
            EFIT.rmaxis = mdsvalue([g_path,'rmaxis']); % magnetic axis [m]
            EFIT.zmaxis = mdsvalue([g_path,'zmaxis']); % magnetic axis [m]
            EFIT.limloc = mdsvalue([a_path,'limloc']); % shape tag
            EFIT.ssep = mdsvalue([a_path,'ssep'])/100; % [m]
            EFIT.beta_N = mdsvalue([a_path,'betan']);
            EFIT.beta_p = mdsvalue([a_path,'betap']); % beta poloidal
            EFIT.kappa = mdsvalue([a_path,'eout']); % elongation
            EFIT.li = mdsvalue([a_path,'li']); % internal inductance
            EFIT.upper_gap = mdsvalue([a_path,'otop'])/100; % [m]
            EFIT.lower_gap = mdsvalue([a_path,'obott'])/100; % [m]
            EFIT.q0 = mdsvalue([a_path,'q0']); % safety factor on axis
            EFIT.qstar = mdsvalue([a_path,'qstar']);
            EFIT.q95 = mdsvalue([a_path,'q95']); % safety factor @ 95% flux
            EFIT.V_loop = mdsvalue([a_path,'vloopt']); % loop voltage [V]
            EFIT.Wmhd = mdsvalue([a_path,'wplasm']); % stored energy [J]
            EFIT.lim = mdsvalue([g_path,'lim']); % limiter geometry
            EFIT.bdry = mdsvalue([g_path,'bdry']); % plasma boundary (psin=1)
            mdsclose();
            loaded = true;
            break
        else;
            disp(['ERROR: No psirz data stored in ', efit_trees{i}, ...
                ' for pulse #',num2str(pulse)])
            mdsclose();
        end
    else
        disp(['Error opening ', efit_trees{i}, ...
            ' tree for pulse #',num2str(pulse)])
    end
end
if ~loaded; EFIT = 0; return; end % return no output from the function if no data
mdsdisconnect(); % Remove MDS+ connection

%% Ensure time is on first dimension for retrieved variables

if size(EFIT.time,1)<size(EFIT.time,2); EFIT.time = EFIT.time'; end
if size(ssimag,2)==length(EFIT.time); ssimag = ssimag'; end
if size(ssibry,2)==length(EFIT.time); ssibry = ssibry'; end
if size(EFIT.rhovn,2)==length(EFIT.time); EFIT.rhovn = EFIT.rhovn'; end
if size(EFIT.aminor,2)==length(EFIT.time); EFIT.aminor = EFIT.aminor'; end
if size(EFIT.rmaxis,2)==length(EFIT.time); EFIT.rmaxis = EFIT.rmaxis'; end
n = find(size(psirz)==length(EFIT.time)); % find time dimension for psirz
if length(n)>1; n = 3; end % if more than one, assume time is 3rd dimension
psirz = shiftdim(psirz,n-1); % put time dimension in front (if not already)

%% Normalize the poloidal flux grid (0=magnetic axis, 1=boundary)
% [Translated from D. Eldon's OMFITeqdsk read_basic_eq_from_mds() function] 

psi_norm_f = ssibry - ssimag;
problems = find(psi_norm_f==0);
% Prevent divide by 0 error by replacing 0s in the denominator
psi_norm_f(problems) = 1;
EFIT.psin = (psirz - ssimag)./psi_norm_f;
EFIT.psin(problems,:,:) = 0;
maxpsin = max(max(EFIT.psin,[],2),[],3); % Max value at any point in space for each time slice
filtered = find(maxpsin ~= 0);
EFIT.psin = EFIT.psin(filtered,:,:);
EFIT.time = EFIT.time(filtered); % Filtered EFIT doesn't have blank (all zero) slices.
EFIT.rmaxis = EFIT.rmaxis(filtered);
EFIT.zmaxis = EFIT.zmaxis(filtered);
EFIT.rhovn = EFIT.rhovn(filtered,:);
end