function shot = argon_data(first_run,last_run,T) %time)
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

% For each shot in shotlist, search for intensity pulse of at least
% T seconds in duration and 1000 counts/s in intensity. This 
% constitutes a 'good' Argon gas concentration. Set all elements in 
% shotlist that do not satisfy this criteria to 0, and then eliminate
% them from the list. 

%T = 0.2;
h = waitbar(0,'Step 1/2: Please wait...');
node = '.hirex_sr.analysis.a';
for i=1:n
    [~, status] = mdsopen('spectroscopy',shotlist(i));
    if (mod(status,2) == 1); %#ok<*ALIGN>
        intensity = mdsvalue([node ':int']);
        [~, status] = mdsvalue(['dim_of(' node ':int)']);
        if (mod(status,2) == 1);
            [~, status] = mdsvalue(['dim_of(' node ':vel)']);
            if (mod(status,2) == 1);
     %----------------------------------------------------------%
                indices = find(intensity>1000 & intensity <10000);
                goodtime = length(indices)*.02;
                if goodtime < T        % Note: Edit this section
                    shotlist(i)=0;     % later if you find a better
                end                    % way to identify 'good'
                                       % Argon puffs. 
     %----------------------------------------------------------%
            end
        else
            shotlist(i)=0;
        end
        mdsclose;
    else
        shotlist(i)=0;
    end
    waitbar(i/n);
end

shotlist = shotlist(shotlist~=0);
shot = size(shotlist);

% Declare "shot" structure array with N elements, each representing
% a shot in the shotlist variable above with a 'good' Argon puff.
% Each shot has 6 fields as follows:
%
% shot(i).i_time = 1xn time domain array corresponding to intensity
% shot(i).intensity = 1xn array with Argon count rate (cnts/sec)
% shot(i).velocity = 1xm array with Argon toroidal velocity values
% shot(i).number = shot number
% shot(i).status = 1 if shot done on locked mode calibration day or
%                  0 if otherwise
%
% All arrays are initialized to NaN and string is initially blank.

%m=max(size(times_for_db));

shot(size(shotlist,1))=struct('i_time',NaN(1,125),'intensity', ...
    NaN(1,125), 'velocity', NaN(1,125), ...
    'number',NaN,'status',NaN);

list = importdata('lock_mode_calib_shots.txt');

% Assign intensity data and shot numbers to each shot structure

waitbar(0,h,'Step 2/2: Please wait...');
n = size(shotlist,1);

for i=1:n
    mdsopen('spectroscopy',shotlist(i));
    [v_time, ~] = mdsvalue(['dim_of(' node ':vel)']);
    [i_time, ~] = mdsvalue(['dim_of(' node ':int)']);
    vel = mdsvalue([node ':vel']);
    velocity = interp1(v_time, vel, time, 'linear');
    intensity = mdsvalue([node ':int']);
    intensity(abs(intensity)>10000)=NaN; % Exclude intensity outliers
    velocity(abs(velocity)>200) = NaN; % Exclude velocity outliers
    shot(i).i_time = i_time;
    shot(i).intensity = intensity;
    shot(i).velocity = velocity * 1000; % Convert from km/s to m/s
    shot(i).number = shotlist(i);
    run = floor(shotlist(i)/1000);
    shot(i).status = any(list == run);
    waitbar(i/n);
end
%}
close(h)
