function [flattop_start, flattop_stop] = get_flattop_times(shot);

% Get the programmed Ip signal from the PCS_EAST tree

mdsopen('pcs_east', double(shot));
[ip_prog, status] = mdsvalue('\lmtipref');
ip_prog = ip_prog * 1.e6; % convert from MA to A
ip_prog_time = mdsvalue('dim_of(\lmtipref)');

% For shots before year 2014, the LMTIPREF timebase needs to be shifted by
% 17.0 ms 

if (shot < 44432);
  ip_prog_time = ip_prog_time - 0.0170;
end;

% Now determine the start and stop times of the flattop

dipprog_dt = gradient(ip_prog, ip_prog_time);

flattop_indices = find(abs(dipprog_dt) <= 1.e3); % 1 kA/s limit
flattop_start = ip_prog_time(min(flattop_indices));
flattop_stop  = ip_prog_time(max(flattop_indices));
