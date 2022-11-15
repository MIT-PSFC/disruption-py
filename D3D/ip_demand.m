db = set_database('d3drdb');

ps = db.prepareStatement('select count(1) as "nshots" from disruptions');
rs = ps.executeQuery();
rs.next();
nshots = rs.getObject(1);

shots = zeros(1, nshots);
tdis = zeros(1, nshots);
ip = zeros(1, nshots);
ip_pcs = zeros(1, nshots);

ps = db.prepareStatement(['select shot,ip,t_disrupt ' ...
  'from disruptions order by shot']);
rs = ps.executeQuery();

for ishot = 1:nshots;
  rs.next();
  shots(ishot) = rs.getObject(1);
  ip(ishot) = rs.getObject(2);
  tdis(ishot) = rs.getObject(3);
end;

db.close();
ps.close();
rs.close();
clear db ps rs;

mdsconnect('atlas.gat.com');

for ishot = 1:nshots;
  shot = shots(ishot);
  fprintf(1,'Processing shot %6i (%5i/%5i  %6.2f%%)\n', shot, ishot, ...
    nshots, ishot/nshots*100);
  [iptipp, status] = mdsvalue(['ptdata("iptipp",' num2str(shot) ')']);
  [time  , status] = mdsvalue(['dim_of(ptdata("iptipp",' ...
                               num2str(shot) '))']);
  if (mod(status,2) == 1);
    [dt, indx] = min(abs(time/1.e3 - (tdis(ishot)-0.050)));
    if (dt <= 0.020);
      ip_pcs(ishot) = iptipp(indx);
    else;
      ip_pcs(ishot) = NaN;
    end;
  else;
    ip_pcs(ishot) = NaN;
  end;
end;

mdsdisconnect;
clear iptipp status time dt indx shot ishot ans;

%save('ip_pcs_data.mat', 'shots', 'tdis', 'ip', 'ip_pcs');
