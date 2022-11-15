load('ip_pcs_data.mat'); % This loads: ip, ip_pcs, shots, tdis

% The plasma current can be positive or negative, but the pcs programming
% is only positive.  (The buswork has to be switched around to change the
% polarity of the plasma current.)  Therefore the absolute value of Ip has
% to be used when determining the difference from the programmed demand.

ip_deficit = ip_pcs - abs(ip);

db = set_database('d3drdb');

for i = 1:numel(shots);
    sql_string = ['update disruptions set ip_deficit = ? ' ...
                  ' where shot = ' ? ];
    ps = db.prepareStatement(sql_string);
    ps.setObject(1, ip_deficit(i));
    ps.setObject(2, shots(i));
    rs = ps.executeUpdate();
end;

db.close();
ps.close();
rs.close();
clear db ps rs;
