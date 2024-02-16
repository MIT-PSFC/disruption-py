function data = extract_data(var)

%------------------------------------------------------------------------%
% Creates cell array 'data' with m rows (where m = the number of shots 
% in the disruption warning database with non-empty velocity data) and 
% 5 columns - shot, time, time_until_disrupt, var, string with name of
% var (in that order).
%
% Author: Kevin Montes
% Date Created: 04/14/2017
%
%------------------------------------------------------------------------%

% Add path to Bob Granetz's folder for matlab database commands 
addpath('/home/granetz/matlab')
addpath('/home/granetz/JRT_2016/disruption_warning_database')

% Extract all relevant data from the logbook into cell array 'result'
db = set_database('logbook');
get_columns(db,'disruption_warning');
result = fetch(db,['select shot, time, time_until_disrupt, ',var, ...
    ' from disruption_warning']);

% Extract data from 'result' and put in usable format
shots = int64(cell2mat(result(:,1))); 
t = cell2mat(result(:,2));
t_until_disrupt = cell2mat(result(:,3));
x = cell2mat(result(:,4));
ushots = unique(shots);
n = size(ushots,1);

% Not all shots will have velocity data (the v_0 vector may only 
% have NaN values) since the get_rotation_velocity.m routine filtered 
% for shots that were calibrated and had a minimum argon intensity pulse.
% Therefore, I only want to store shots in the 'data' cell array that have
% velocity data to plot (for 2015 campaign, ~150 shots). 
%
% First, declare empty cell array C and track all 'ushots' indices that 
% have empty velocity vectors. Populate only rows of C that correspond to
% shots with non-empty velocity vectors. Then, equate the 'data' cell array
% to the relevant part of C. 

C = cell(n,5);
C(:,5) = {var};
indices = 1:n;

for i=1:n
    bool = ismember(shots,ushots(i));
    if all(isnan(x(bool)))
        indices(i) = 0;
    else
        C(i,1) = {ushots(i)};
        C(i,2) = {t(bool)};
        C(i,3) = {t_until_disrupt(bool)};
        C(i,4) = {x(bool)};
    end
end

data = C(indices(indices~=0),:);