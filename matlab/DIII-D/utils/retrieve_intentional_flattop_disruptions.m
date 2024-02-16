db = set_database('d3drdb');
result = fetch(db, ['select shot, time, time_until_disrupt, ip,' ...
  'ip_error, dipprog_dt, power_supply_railed,intentional_disruption,' ...
  'other_hardware_failure from disruption_warning order by shot,time']);
close(db);

shot = int32(cell2mat(result(:,1)));
time = cell2mat(result(:,2));
time_until_disrupt = cell2mat(result(:,3));
ip = cell2mat(result(:,4));
ip_error = cell2mat(result(:,5));
dipprog_dt = cell2mat(result(:,6));
power_supply_railed = cell2mat(result(:,7));
intentional_disruption = cell2mat(result(:,8));
other_hardware_failure = cell2mat(result(:,9));

clearvars -except shot time time_until_disrupt ip ...
  ip_error dipprog_dt power_supply_railed intentional_disruption ...
  other_hardware_failure;

ipprog = ip - ip_error;
D
indices_disrupt = find(~isnan(time_until_disrupt));
indices_intentional_flattop = find(abs(dipprog_dt) <= 2.e3 & ...  % 2 kA/s limit
  abs(ip_prog) > 100e3 &  intentional_disruption == 1);
indices_disrupt_time = find(time_until_disrupt == 0);

indices_intentional_flattop_disruptions = find(abs(dipprog_dt) <= 2.e3 & ...
  abs(ip_prog) > 100e3 &  intentional_disruption == 1 & ...
  ~isnan(time_until_disrupt));

intentional_disruptions = unique(shot(indices_intentional_flattop_disruptions));
% 325

duration = zeros(length(intentional_disruptions),1);
intentional_flattop_time = time(indices_intentional_flattop_disruptions);

for ishot = 1:length(intentional_disruptions);
	fprintf(1,'shot %5d\n',intentional_disruptions(ishot));
	idx = find(shot(indices_intentional_flattop_disruptions)== intentional_disruptions(ishot));
	tt = intentional_flattop_time(idx);
	duration(ishot) = tt(end)-tt(1);	
end;

long = find(duration > 1); %supposedly enough to get into Hmode?
shotlist =  intentional_disruptions(long);

