%% Grab shotlist and DB info
addpath('/home/granetz/JRT_2016/disruption_warning_database')
load('/home/montes/dataset_partition.mat')
shotlist = [testShotlist;trainingShotlist;validationShotlist];
db = set_database('logbook');
result = fetch(db,['select shot, time, ip, ip_error from disruption_warning ' ...
    'order by shot, time'],'DataReturnFormat','cellarray');
shot = int32(cell2mat(result(:,1)));
time = cell2mat(result(:,2));
ip = cell2mat(result(:,3));
ip_error = cell2mat(result(:,4));

%% Loop through shots
N = length(shotlist);
ip_prog = cell(length(shotlist),2);
for i=1:N
    shot_num = shotlist(i);
    indices = find(shot==shot_num);
    ip_prog{i,1} = time(indices);
    ip_prog{i,2} = ip(indices)-ip_error(indices);
    if ~mod(i,100)
        disp([num2str(i) '/' num2str(N)])
    end
end