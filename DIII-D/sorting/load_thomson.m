function TS_out = load_thomson(pulse,sig_type,varargin)
%LOAD_THOMSON Retrieves DIII-D Thomson scattering data from DIII-D
%   Author: Kevin Montes
%   Date: March 2019

mdsconnect('atlas.gat.com');

omit_last = false; % if true, omit the last Thomson channel 

if ~strcmp(sig_type,'temp') & ~strcmp(sig_type,'density')
    disp(['ERROR: "',sig_type,'" not a valid signal type for Thomson'])
    return
end

% Get names of laser systems to use
if length(varargin) > 0
    TS_systems = varargin{1};
else
    TS_systems = {'core','tangential'};
end

% Choose whether data used is 'raw' or 'blessed' by Thomson group
version = 'ptdata'; % set to blessaed, unblessed, or ptdata
if strcmp(version,'blessed')
    mds_path = '\top.ts.blessed.';
elseif strcmp(version,'unblessed')
    mds_path = '\top.ts.revisions.revision00.';
elseif strcmp(version,'ptdata')
	mds_path = '\top.ts.blessed.';
	sig = 'ne'; if strcmp(sig_type,'temp'); sig='te'; end
	for i=1:length(TS_systems); suffix.(TS_systems{i}) = TS_systems{i}(1:3); end
	c0 = 0; cf = 1; % start channel indexing at 1
	format = '%02d'; % add 1 trailing 0 on channel number
end

% Account for change in pointname formatting in 2017... Use old format if old shot
if pulse<172749 % first shot on Sep. 19, 2017
	suffix.tangential = 'hor'; % 'horizontal' suffix used in old naming system
	c0 = 1; cf = 0;
	format = '%01d'; % no trailing 0 on channel number
end

% Open electrons tree and retrieve all relevant data from each subsystem
results = cell(1,2);
[shotopened, status]=mdsopen('electrons', double(pulse));
if (mod(status,2)==1);
    for i = 1:length(TS_systems)
        laser = TS_systems{i};
		subtree = [mds_path,laser];
    	[time, status] = mdsvalue(['dim_of(',subtree,':',sig_type,',0)']);
	    if (mod(status,2)==1)
	        TS.R = mdsvalue([subtree,':r']);
	        TS.z = mdsvalue([subtree,':z']);
			if strcmp(version,'ptdata')
				TS.time = mdsvalue(['dim_of(ptdata("tss',suffix.(laser),sig, ...
					sprintf(format,c0),'"))'])/1e3;
				TS.data = NaN(length(TS.time),length(TS.R));
				channels = [c0:length(TS.R)-cf];
				for j=1:length(channels)
					TS.data(:,j) = mdsvalue(['ptdata("tss',suffix.(laser), ...
										sig,sprintf(format,channels(j)),'")']);
				end
			else
	        	TS.data = mdsvalue([subtree,':',sig_type]);
				TS.time = time/1.e3; % convert to ms to s
			end
	        if size(TS.R,2)==1; TS.R = TS.R'; end
	        if size(TS.z,2)==1; TS.z = TS.z'; end
	        if size(TS.time,1)==1; TS.time = TS.time'; end
	        if size(TS.data,1)~=length(TS.time); TS.data = TS.data'; end
	        % Place NaNs for broken channels
	        TS.data(find(TS.data==0)) = NaN;
	        results{i} = TS;
	    else
	        disp(['Error loading TS ',laser, ...
	            ' data for shot ',num2str(pulse)])
	    end
    end
    mdsclose();
else;
    disp(['Error opening electrons tree!'])
end

% If both systems/lasers available, combine them and interpolate the data
% from the tangential system onto the finer (core) timebase
if ~isempty(results{1})
    if ~isempty(results{2}) % both core & tangential systems are available
        TS_out.data = [results{1}.data, ...
            interp1(results{2}.time,results{2}.data,results{1}.time)];
        TS_out.R = [results{1}.R,results{2}.R];
        TS_out.z = [results{1}.z,results{2}.z];
        TS_out.time = results{1}.time;
    else % only core system is available
        TS_out = results{1};
    end
    TS_out.mds_path = mds_path;
else
    return % Neither system is available -> return with no output
end

if omit_last
	TS_out.data(:,end) = []; % omits tan09 for consistency, which does not get returned in realtime
	TS_out.R(end) = [];
	TS_out.z(end) = [];
end

if all(all(isnan(TS_out.data)))
    disp('All Thomson channels empty... Aborting!');
    clear
    return
end

mdsdisconnect();

end
