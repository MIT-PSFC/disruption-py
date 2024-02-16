db = set_database('d3drdb'); % My routine to connect to database

ps = db.prepareStatement(['alter table disruptions add ip_programmed ' ...
                          'float null']);
rs = ps.executeUpdate();

db.close();
ps.close();
rs.close();
clear db ps rs
