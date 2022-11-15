function TS_out = get_thomson(pulse,varargin)
%LOAD_THOMSON Retrieves DIII-D Thomson scattering data from DIII-D
%   Author: Kevin Montes
%   Date: March 2019

mdsconnect('atlas.gat.com');

if nargin > 1
    TS_systems = varargin{1};
else
    TS_systems = {'core','tangential'};
end

results = cell(1,2);
[shotopened, status]=mdsopen('electrons', double(pulse));
if (mod(status,2)==1);
    for i = 1:length(TS_systems)
        laser = TS_systems{i};
        subtree = ['\top.ts.blessed.',laser,':'];
        [time, status] = mdsvalue(['dim_of(',subtree,'temp,0)']);
        if (mod(status,2)==1)
            TS.Te = mdsvalue([subtree,'temp']);
	    TS.ne = mdsvalue([subtree,'density']);
            TS.R = mdsvalue([subtree,'r']);
            TS.z = mdsvalue([subtree,'z']);
            TS.time = time/1.e3; % convert to ms to s
            if size(TS.R,2)==1; TS.R = TS.R'; end
            if size(TS.z,2)==1; TS.z = TS.z'; end
            if size(TS.time,1)==1; TS.time = TS.time'; end
            if size(TS.Te,1)~=length(TS.time); TS.Te = TS.Te'; end
	    if size(TS.ne,1)~=length(TS.time); TS.ne = TS.ne'; end
            % Place NaNs for broken channels
	    TS.Te(find(TS.Te==0)) = NaN;
            TS.ne(find(TS.ne==0)) = NaN; 
            if strcmp(laser,'core')
                % Adjust for broken core channel 41
                TS.Te = TS.Te(:,1:end-1);
		TS.ne = TS.ne(:,1:end-1);
                TS.R = TS.R(1:end-1);
                TS.z = TS.z(1:end-1);
            end 
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
        TS_out.Te = [results{1}.Te, ...
            interp1(results{2}.time,results{2}.Te,results{1}.time)];
	TS_out.ne = [results{1}.ne, ...
            interp1(results{2}.time,results{2}.ne,results{1}.time)];
        TS_out.R = [results{1}.R,results{2}.R];
        TS_out.z = [results{1}.z,results{2}.z];
        TS_out.time = results{1}.time;
    else % only core system is available
        TS_out = results{1};
    end
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
