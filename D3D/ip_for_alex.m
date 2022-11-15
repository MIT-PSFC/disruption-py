function ip_data = ip_for_alex(shotlist);

mdsconnect('atlas.gat.com');
ip_data = [];

% Start for loop over the input shots, acquire ip signal at each
for i=1:length(shotlist)
    [shotopened, status]=mdsopen('d3d', shotlist(i));
    ip = mdsvalue(['ptdata("ip", ' num2str(shotlist(i)) ')']);
    [iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' num2str(shotlist(i)) '))']);
    if (mod(status,2) == 1);
        iptime = iptime/1.e3; % convert ms to s
        if size(ip,2)>1; ip = ip'; end
        if size(iptime,2)>1; iptime = iptime'; end
        shot = ones(size(iptime))*shotlist(i);
        shot_data = [ip,iptime,shot];
        ip_data = [ip_data;shot_data];
    end
end

end

        
        
        
    
