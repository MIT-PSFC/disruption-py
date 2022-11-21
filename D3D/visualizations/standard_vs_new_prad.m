warning_status = warning;
warning('off')

db_retrieve = set_database('d3drdb');
variables = {'time'; 'time_until_disrupt'; 'shot'};
max_charlength = size(char(variables),2);

for ivars = 1:length(variables);
  varname = char(variables(ivars));
  charlength = length(varname);
  strng = [varname blanks(max_charlength - charlength)];
  fprintf(1, 'retrieving %s  (%2i/%2i)\n', strng, ivars, length(variables));
  result_retrieve = fetch(db_retrieve, ['select ' varname ...
    ' from disruption_warning order by shot, time']);
  eval([varname '= cell2mat(result_retrieve);']);
end;

shot_db = int32(shot); % convert shot numbers to integers
close(db_retrieve);
% Return warning message status to original setting
if (strcmpi(warning_status(1).state, 'on'));
  warning('on');
end;

clearvars db_retrieve varname ivars result_retrieve max_charlength ...
          charlength strng warning_status variables shot;

load('flattop_disruptions2015.mat');
%shot = 164716;
mdsconnect('atlas.gat.com');

nshots = length(shotlist); 
ndigits = int2str(floor(log10(nshots)) + 1); % (for formatting in the 
                                             %  upcoming fprintf statement)
smoothing_window = 0.025; % use 25 ms causal smoothing window

sum_prad_new = NaN(length(shotlist), 1);
sum_prad_tot = NaN(length(shotlist), 1);

for ishot=1:nshots;
  shot = shotlist(ishot);
  index = find(shot_db == shot);
  tmp = time_until_disrupt(index);
  time_tmp = time(index);
  
  
  fprintf(1,['Processing shot %7i  (%' ndigits 'i/%' ndigits ...
    'i  %6.2f%%)\n'], shot, ishot, nshots, ishot/nshots*100);

    mdsopen('bolom', shot);
    prad_tot = mdsvalue('\prad_tot');
    prad_tot_timebase = mdsvalue('dim_of(\prad_tot)')/1.e3;
    mdsclose;

    a_structure = getbolo(shot, smoothing_window*1.e3);
    b_structure = powers(a_structure);
    prad_new = b_structure.pwrmix; % watts
    prad_new_timebase = a_structure.rawtime; % seconds
    
    mdsopen('d3d', shot);
    ip = mdsvalue(['ptdata("ipprobesf", ' num2str(shot) ')']);
    [iptime, status] = mdsvalue(['dim_of(ptdata("ipprobesf", ' num2str(shot) '))']);
    if (mod(status,2) == 1);
        iptime = iptime/1.e3; % convert ms to s
    end;
    ipprog = mdsvalue(['ptdata("iptipp", ' num2str(shot) ')']);
    [ipprog_time, status] = mdsvalue(['dim_of(ptdata("iptipp", ' num2str(shot) '))']);
    if (mod(status,2) == 1);
        ipprog_time = ipprog_time/1.e3; % convert ms to s
    end;
    mdsclose;
    
    ipprog = interp1(ipprog_time, ipprog, iptime, 'linear');
    
    time_end = time_tmp(tmp ==0);
    time_end_std = time_end+20.e-3;
    %find average value of ip in flattop
    %d_ip = GSASTD(iptime(indices_flattop),ip(indices_flattop),1,20,3,1,0);
    
    indices_flattop = find(abs(gradient(ipprog,iptime)) <= 2.e3 & ...
                    abs(ipprog) > 100e3 & iptime < time_end); % 2 kA/s limit & power supplies not railed

    %dip_time = iptime(indices_flattop);
    %time_start = time_end - 10.e-3%max(iptime(find(ip==max(ip(indices_flattop)))))
    %time_start2 =  dip_time(find(d_ip==max(d_ip)))
    time_start = iptime(ip==max(ip(indices_flattop)));
    if (length(time_start)>1) | (diff(time_start)<2.e-3)
        time_start = mean(time_start);
    end
    time_start_std = time_start -20.e-3;
    
    index_prad_new = find(prad_new_timebase >= time_start & ...
                    prad_new_timebase <= time_end);
    index_prad_tot = find(prad_tot_timebase >= time_start_std & ...
                    prad_tot_timebase <= time_end_std);

    sum_prad_new(ishot) = sum(prad_new(index_prad_new))/1.e9;
    sum_prad_tot(ishot) = sum(prad_tot(index_prad_tot))/1.e9;

    %{
    figure()
    plot(iptime,ip/1.e6,'b-')
    hold on
    plot(iptime(find(iptime>=time_start & iptime <= time_end)),...
        ip(find(iptime>=time_start & iptime <= time_end))/1.e6,'g-')
    hold on
    plot(iptime(find(iptime>=time_start2 & iptime <= time_end)),...
        ip(find(iptime>=time_start2 & iptime <= time_end))/1.e6,'r-')
    hold on
    plot(iptime(indices_flattop), ipprog(indices_flattop)/1.e6,'k-')
    %}
end;
mdsclose;

%figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
%       'PaperPositionMode', 'auto');
plot(sum_prad_new, sum_prad_tot ,'.b', 'markersize', 13);
hold on;
plot([0 xlim], [0 xlim], '-r', 'linewidth', 1.5);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
xlabel('sum(P_{rad} new) [GW]', 'FontSize', 14);
ylabel('sum(P_{rad} std) [GW]', 'FontSize', 14);

%print('my_standard_prad_vs_new_prad.png', '-dpng');
