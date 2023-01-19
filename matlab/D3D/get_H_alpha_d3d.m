function H_alpha=get_H_alpha_d3d(shot, timebase)
if (size(timebase,2) > 1)
  timebase = transpose(timebase);
end

H_alpha = NaN(size(timebase));

mdsconnect('atlas.gat.com');
[~, status] = mdsopen('d3d', shot);
if (mod(status,2) == 0)
    mdsclose;
    return
end

[H_alpha_time, status]=mdsvalue('dim_of(\fs04)');
if (mod(status,2) == 1 && length(H_alpha_time) > 1)
    H_alpha_time = H_alpha_time/1.e3;
    H_alpha = mdsvalue('\fs04');
    H_alpha = interp1(H_alpha_time, H_alpha, timebase, 'linear');
else
    mdsclose;
    return;
end
mdsclose;
end