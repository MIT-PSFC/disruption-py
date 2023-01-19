shotlist = NaN(1,250000);

db = set_database('d3drdb');
ps = db.prepareStatement('select shot from disruptions order by shot');
rs = ps.executeQuery();
counter = 0;
while (rs.next());
  counter = counter + 1;
  disruption_shot = rs.getObject(1);
  shotlist(counter) = disruption_shot;
end;
db.close();
ps.close();
rs.close();
