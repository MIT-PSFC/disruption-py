function EFIT = load_efit(pulse,varargin)
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

mdsconnect('atlas.gat.com'); % Start MDS+ connection for DIII-D

% List EFIT trees to try, in order of priority
if length(varargin) > 0
    efit_trees = varargin{1};
else
    efit_trees = {'efit01','efit02','efit03','efit04'};
end

%% Pull the relevant data from MDS+
g_path = '\top.results.geqdsk:';
a_path = '\top.results.aeqdsk:';
i = 0; loaded = false;
while i < length(efit_trees)
    i = i+1;
    % Open EFIT tree for this shot and grab the grid
    [shotopened, status]=mdsopen(efit_trees{i}, double(pulse));
    if mod(status,2)==1;
        [EFIT.time,status] = mdsvalue(['dim_of(',g_path,'psirz,2)']);
        if mod(status,2)==1;
	    	EFIT.tree = efit_trees{i};
            EFIT.time = EFIT.time/1e3; % [s]
            EFIT.z = mdsvalue([g_path,'z']); % [m]
            EFIT.r = mdsvalue([g_path,'r']); % [m]
            rhovn = mdsvalue([g_path,'rhovn']); % normalized radius
            psirz = mdsvalue([g_path,'psirz']); % measured poloidal flux
			qpsi = mdsvalue([g_path,'qpsi']); % safety factor on uniform flux grid
            ssimag = mdsvalue([g_path,'ssimag']);
            ssibry = mdsvalue([g_path,'ssibry']);
            EFIT.aminor = mdsvalue([a_path,'aminor']); % minor radius [m]
            EFIT.aminor_time = mdsvalue(['dim_of(',a_path,'aminor)'])/1e3;
            EFIT.rmaxis = mdsvalue([g_path,'rmaxis']);
            EFIT.zmaxis = mdsvalue([g_path,'zmaxis']);
            EFIT.limloc = mdsvalue([a_path,'limloc']);
			EFIT.Wmhd = mdsvalue([a_path,'wmhd']);
			EFIT.Wdia = mdsvalue([a_path,'wdia']);
            EFIT.ssep = mdsvalue([a_path,'ssep']);
			EFIT.betap = mdsvalue([a_path,'betap']);
			EFIT.li = mdsvalue([a_path,'li']);
			EFIT.q95 = mdsvalue([a_path,'q95']);
			EFIT.H98 = mdsvalue(['\top.results.confinement.times.scalings:h_thh98y2']);
			EFIT.kappa = mdsvalue([a_path,'kappa']);
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
if size(ssimag,2)==length(EFIT.time); ssimag = ssimag'; end
if size(ssibry,2)==length(EFIT.time); ssibry = ssibry'; end
if size(rhovn,2)==length(EFIT.time); rhovn = rhovn'; end
if size(qpsi,2)==length(EFIT.time); qpsi = qpsi'; end
n = find(size(psirz)==length(EFIT.time)); % find time dimension for psirz
if length(n)>1; n = 3; end % if more than one, assume time is 3rd dimension
psirz = shiftdim(psirz,n-1); % put time dimension in front (if not already)

%% Normalize the poloidal flux grid (0=magnetic axis, 1=boundary)
% [Translated from D. Eldon's OMFITeqdsk read_basic_eq_from_mds() function] 
psi_norm_f = ssibry - ssimag;
problems = find(psi_norm_f==0);
psi_norm_f(problems) = 1; % Prevent divide by 0 error by replacing 0s in the denominator
psin = (psirz - ssimag)./psi_norm_f;
psin(problems,:,:) = 0;
maxpsin = max(max(psin,[],2),[],3); % Max value at any point in space for each time slice
filtered = find(maxpsin ~= 0);

%% Format supplemental retrieved variables
names = fieldnames(EFIT);
for j=1:length(names)
	if isrow(EFIT.(names{j})) & length(EFIT.(names{j}))==length(EFIT.time)
 		EFIT.(names{j}) = EFIT.(names{j})'; % make 1st dim time
		if ~ischar(EFIT.(names{j}));
			EFIT.(names{j}) = EFIT.(names{j})(filtered); % filter out bad EFIT slices
		end
	end
end
EFIT.psin = psin(filtered,:,:);
EFIT.rhovn = rhovn(filtered,:);
EFIT.qpsi = qpsi(filtered,:);

end
