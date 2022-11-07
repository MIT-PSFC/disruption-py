db = set_database('logbook');
get_columns(db,'disruption_warning');
result = fetch(db, ['select shot from disruption_warning order by shot']);
shot = int32(cell2mat(result(:,1)));
shotlist = unique(shot);
n = length(shotlist);

% Create 9x4 cell array, where 9 rows contain data for each of the 9 ECE
% channels, and the 4 columns contain the following data:
%
%   Column 1    Column 2    Column 3    Column 4
%      T         T_time       rad        rad_time

Te_data = cell(9,4);
node = '.ece.gpc_results';

h = waitbar(0,'Please wait ...');
shot = 115
    [~, status] = mdsopen('electrons',shotlist(i));
    if (mod(status,2) == 1)
        for j=1:9
            Te_data{j,1} = mdsvalue([node '.te:te' num2str(j)]);
            [Te_data{j,2},status] = mdsvalue(['dim_of(' node '.te:te'...
                num2str(j) ')']);
            Te_data{j,3} = mdsvalue([node '.rad:r' num2str(j)]);
            [Te_data{j,4},status] = mdsvalue(['dim_of(' node ...
                '.rad:r' num2str(i) ')']);
        end
    end
    waitbar(i/n,h,['Shot ' num2str(i) ' of ' num2str(n)]); 
end
close(h)