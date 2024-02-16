function smeared_times = smear(times, dt);

if (nargin > 1);
  delta = dt;
else;
  delta = 0.001;
end;
smeared_times = times;

indices = find(abs( round(times*1.e3)/1.e3 - times ) < 1e-6);
n_smear = length(indices);
rng('default');
smeared_dt = -delta/2 + delta*rand(n_smear,1);
smeared_times(indices) = smeared_times(indices) + smeared_dt;

end