db = set_database('logbook');
result = fetch(db,['select shot,time_until_disrupt,ip,ip_error,' ...
                'dipprog_dt,time from disruption_warning order by ' ...
                'shot,time']);
close(db);

shots = int32(cell2mat(result(:,1)));
tdis = cell2mat(result(:,2));
ip = cell2mat(result(:,3));
iperr = cell2mat(result(:,4));
didt = cell2mat(result(:,5));
times = cell2mat(result(:,6));

% Take account of the Ip polarity.  (Maybe I should add a 'polarity'
% parameter to the database?)

unique_shots = unique(shots);
for i = 1:length(unique_shots);
  indices = find(shots == unique_shots(i));
  polarity = sign(sum(ip(indices)));
  iperr(indices) = iperr(indices) * polarity;
  didt(indices) = didt(indices) * polarity;
end;

flattop  = find(didt == 0);
rampdown = find(didt  < 0);
rampup   = find(didt  > 0);

figure('menubar','none','toolbar','none');
plot(-tdis(flattop)*1.e3, iperr(flattop)/1.e3, '.b');
xlim([-350,0]);
ylim([-100,100]);
xlabel('Time before disruption [ms]', 'fontsize', 12);
ylabel('Ip error ( = Ip - Ip_{programmed}) [kA]', 'fontsize', 12);
title('During flattop', 'fontsize', 14);
hold on;
plot([-350,0],[0,0],'g','linewidth',2);
hold off;

figure('menubar','none','toolbar','none');
plot(-tdis(rampdown)*1.e3, iperr(rampdown)/1.e3, '.b');
xlim([-350,0]);
ylim([-100,100]);
xlabel('Time before disruption [ms]', 'fontsize', 12);
ylabel('Ip error ( = Ip - Ip_{programmed}) [kA]', 'fontsize', 12);
title('During rampdown', 'fontsize', 14);
hold on;
plot([-350,0],[0,0],'g','linewidth',2);
hold off;

figure('menubar','none','toolbar','none');
plot(-tdis(rampup)*1.e3, iperr(rampup)/1.e3, '.b');
xlim([-350,0]);
ylim([-100,100]);
xlabel('Time before disruption [ms]', 'fontsize', 12);
ylabel('Ip error ( = Ip - Ip_{programmed}) [kA]', 'fontsize', 12);
title('During rampup', 'fontsize', 14);
hold on;
plot([-350,0],[0,0],'g','linewidth',2);
hold off;

figure('menubar','none','toolbar','none');
plot(times, iperr/1.e3, '.b');
xlim([0,2.3]);
ylim([-100,100]);
xlabel('Time during shot [s]', 'fontsize', 12);
ylabel('Ip error ( = Ip - Ip_{programmed}) [kA]', 'fontsize', 12);
%title('During rampup', 'fontsize', 14);
hold on;
plot([0,2.3],[0,0],'g','linewidth',2);
hold off;
