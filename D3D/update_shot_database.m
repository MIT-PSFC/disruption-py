% Determine the last shot that is already in the shot database file.  A
% scan of possible new shots should start at the last shot plus one.

read_shot_list;
startshot = shotlist(max(find(~isnan(ipmax) & (~isnan(eos))))) + 1;
endshot = startshot + 25;
shot_database;

read_shot_list;
lastshot = shotlist(max(find(~isnan(ipmax) & (~isnan(eos)))));

% Truncate the records at the end of the file corresponding to shots which
% have not yet occurred.  (i.e. MDS trees have not been created yet, so
% ipmax and EOS are NaN)

sedcmd = sprintf(['sed -e '' ' num2str(lastshot+4) ...
                  ',$d'' -i EAST_shot_list.txt']);
unix(sedcmd,'-echo');

% Read in the up-to-date list of shots, ipmax, eos, timestamps

read_shot_list;
new_lastshot = shotlist(max(find(~isnan(ipmax) & (~isnan(eos)))));

% Now update the disruption database

disruption_database(startshot:new_lastshot);
