% Creates data vector of radiated power after reading in time, shot, and
% p_rad. Uses trapezoidal sum to approximate integral via Matlab's 
% trapz() function. 

shotlist = unique(shot,'stable');
int_p_rad = [];
for i=1:size(shotlist,1)
    indices = find(shot == shotlist(i));
    p = p_rad(indices);
    t = time(indices);
    I = NaN(size(p));
    for j = 2:size(p,1)
        I(j) = trapz(t(1:j),p(1:j));
    end 
    int_p_rad = [int_p_rad;I];
end

clear shotlist indices i j p t I