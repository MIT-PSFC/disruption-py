function sigs = find_signal_path(keyword)

%% This function searches the SIGNAL_NAMES database at DIII-D to find
% the MDSplus paths for signals matching the keyword input argument
%
% Input:
%	- keyword: a string for which to search the signal names in the database
%
% Output:
%	- sigs: an Nx3 cell array with the signal name in the 1st column, the
%			signal type in the 2nd column, and the full signal path in
%			the 3rd column
%
% Author: Kevin Montes
% Date: 02/2020
%
%----------------------------------------------------------------------------

if ~ischar(keyword)
	disp('Input keyword must be a character vector ...')
	return
end

db = set_database('d3drdb');
if ~isopen(db)
	disp('Failed to open SQL database ...')
	return
end

result = fetch(db,['select name, type, full_path from signal_names']);
indx = find(~cellfun(@isempty,strfind(result(:,1),keyword)));
sigs = result(indx,:);

