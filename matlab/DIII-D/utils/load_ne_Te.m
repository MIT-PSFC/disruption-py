function TS_out = load_ne_Te(pulse,data_source,varargin)
%LOAD_THOMSON Retrieves DIII-D Thomson scattering data from DIII-D
%   Author: Kevin Montes
%   Date: March 2019
% 	Inputs:
%		- pulse: shot number to pull data for
%		- version: set to blessed, unblessed, or ptdata
%-----------------------------------------------------------------

mdsconnect('atlas.gat.com');

% Get names of laser systems to use
if nargin > 2
    TS_systems = varargin{1};
else
    TS_systems = {'core','tangential'};
end

% Choose whether data used is 'blessed' by Thomson group (or 'unblessed')
if strcmp(data_source,'blessed')
    mds_path = '\top.ts.blessed.';
elseif strcmp(data_source,'unblessed')
    mds_path = '\top.ts.revisions.revision00.';
elseif strcmp(data_source,'ptdata')
	mds_path = '\top.ts.blessed.';
	for i=1:length(TS_systems); 
		suffix.(TS_systems{i}) = TS_systems{i}(1:3); % get first 3 letters for pointname
	end
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
if (mod(status,2)==1)
    for i = 1:length(TS_systems)
        laser = TS_systems{i};
		subtree = [mds_path,laser];
    	[time, status] = mdsvalue(['dim_of(',subtree,':temp,0)']);
	    if (mod(status,2)==1)
	        TS.R = mdsvalue([subtree,':r']); % major radial position of measurement
	        TS.z = mdsvalue([subtree,':z']); % vertical position of measurement
			if strcmp(data_source,'ptdata')
				TS.time = mdsvalue(['dim_of(ptdata("tss',suffix.(laser),'te', ...
					sprintf(format,c0),'"))'])/1e3;
				TS.Te = NaN(length(TS.time),length(TS.R));
				TS.ne = NaN(length(TS.time),length(TS.R));
				channels = [c0:length(TS.R)-cf];
				for j=1:length(channels)
					TS.Te(:,j) = mdsvalue(['ptdata("tss',suffix.(laser), ...
										'te',sprintf(format,channels(j)),'")']);
					TS.ne(:,j) = mdsvalue(['ptdata("tss',suffix.(laser), ...
										'ne',sprintf(format,channels(j)),'")']);
				end
			else
	        	TS.Te = mdsvalue([subtree,':temp']);
				TS.ne = mdsvalue([subtree,':density']);
				TS.Te_error = mdsvalue([subtree,':temp_e']);
				TS.ne_error = mdsvalue([subtree,':temp_e']);
				TS.time = time/1.e3; % convert to ms to s
			end
			names = fields(TS);
			for j=1:length(names) % make sure all fields have consistent dimensions
				[s1,s2] = size(TS.(names{j}));
				if (s1==length(TS.R)) | (s2==length(TS.time)) 
					TS.(names{j}) = TS.(names{j})'; % transpose to [time x channels] array
				end
			end
	        % Place NaNs for broken channels
	        TS.Te(find(TS.Te==0)) = NaN;
			TS.ne(find(TS.ne==0)) = NaN;
	        results{i} = TS;
	    else
	        disp(['Error loading TS ',laser, ...
	            ' data for shot ',num2str(pulse)])
	    end
    end
    mdsclose();
else
    disp(['Error opening electrons tree!'])
end

% If both systems/lasers available, combine them and interpolate the data
% from the tangential system onto the finer (core) timebase
if ~isempty(results{1})
    if ~isempty(results{2}) % both core & tangential systems are available
		vars = find(~(strcmp(names,'R') | strcmp(names,'z') | strcmp(names,'time')));
		for i=1:length(vars)
			var = names{vars(i)};
		    TS_out.(var) = [results{1}.(var), ...
		        interp1(results{2}.time,results{2}.(var),results{1}.time)];
		end
        TS_out.R = [results{1}.R,results{2}.R];
        TS_out.z = [results{1}.z,results{2}.z];
        TS_out.time = results{1}.time;
    else % only core system is available
        TS_out = results{1};
    end
    TS_out.data_source = data_source;
	TS_out.shot = pulse;
else
    return % Neither system is available -> return with no output
end

if all(all(isnan(TS_out.Te)))
    disp('All Thomson channels empty... Aborting!');
    clear
    return
end

mdsdisconnect();

end
