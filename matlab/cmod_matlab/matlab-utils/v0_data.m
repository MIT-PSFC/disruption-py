function [count,calib_count] = v0_data(first_run,last_run,min) %time)
%{
This function goes through all shots in the run range indicated by
the two input arguments (first_run, last_run) and lists only the
shots for which is there is a significant amount of Argon. 
%}

% First, search the logbook for shots in the range

addpath('/home/granetz/matlab');
db = set_database('logbook');
results = fetch(db,['select distinct(shot) from entries where ' ...
    'run between ' int2str(first_run) ' and ' int2str(last_run)]);
shotlist = int32(cell2mat(results));

% Filter out test shots and store integer shot numbers in shotlist
shotlist2 = shotlist - floor(shotlist/1000)*1000;
shotlist = shotlist(shotlist2 > 0 & shotlist2 < 100);
n = size(shotlist,1);
shotlist = double(shotlist);
shotlist = [shotlist,NaN(size(shotlist))];

% For each shot in shotlist, search for intensity pulse of at least
% T seconds in duration and 1000 counts/s in intensity. This 
% constitutes a 'good' Argon gas concentration. Set all elements in 
% shotlist that do not satisfy this criteria to 0, and then eliminate
% them from the list. 

%T = 0.2;
h = waitbar(0,'Please wait...');
node = '.hirex_sr.analysis.a';

for i=1:n
    [~, status] = mdsopen('spectroscopy',shotlist(i,1));
    if (mod(status,2) == 1); %#ok<*ALIGN>
        intensity = mdsvalue([node ':int']);
        [~, status] = mdsvalue(['dim_of(' node ':int)']);
        if (mod(status,2) == 1);
            % Use 'units_of()' to get units if needed
            [~, status] = mdsvalue(['dim_of(' node ':vel)']);
            if (mod(status,2) == 1);
     %----------------------------------------------------------%
                indices = find(intensity>min & intensity <10000);
                goodtime = length(indices)*.02;
                shotlist(i,2) = goodtime; 
     %----------------------------------------------------------%
            end
        end
        mdsclose;
    end
    waitbar(i/n);
end

% Initialize time and count variables for all and calibrated shots
time = [0:0.01:1.5];
count = NaN(size(time));
calib_count = count;

% Calculate count for all shots
for i=1:size(time,2)
    count(i) = size(find(shotlist(:,2)>time(i)),1);
end

% Calculate count of calibrated shots
list = importdata('lock_mode_calib_shots.txt');
shots = shotlist(:,1);
goodtimes = shotlist(:,2);
flags = zeros(size(shots)); % to flag calibrated shots
for i=1:size(shots)
    flags(i) = any(floor(shots(i)/1000) == list); 
end
goodtimes = goodtimes(find(flags));  %stores only calibrated goodtimes
for i=1:size(time,2)
    calib_count(i) = size(find(goodtimes>time(i)),1);
end

% Plot graph of threshold pulse time vs. number of good shots
plot(time,count,time,calib_count)
legend('All Shots','Calibrated Shots')
xlabel('Threshold Pulse Time (s)')
ylabel('Number of Good Shots')
title('Shots w/ Intensity Pulse Above Threshold')

close(h);