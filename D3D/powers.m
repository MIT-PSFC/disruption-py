function b = powers(a);

kappa=[1.976e8,2.060e8,2.146e8,2.319e8,2.277e8,2.773e8, ...
       2.845e8,2.877e8,2.780e8,2.692e8,2.519e8,2.509e8, ...
       1.004e8,0.935e8,0.893e8,0.790e8,0.711e8,0.683e8, ...
       0.861e8,0.834e8,0.812e8,0.860e8,0.941e8,0.991e8, ...
       1.605e8,1.588e8,1.592e8,1.561e8,1.496e8,0.700e8, ...
       0.698e8,0.699e8,0.719e8,0.774e8,0.865e8,0.692e8, ...
       0.683e8,0.698e8,0.739e8,1.742e8,1.771e8,1.898e8, ...
       1.935e8,1.994e8,2.035e8,2.033e8,2.030e8,0.848e8, ...
       4.632e8,4.370e8,4.230e8,3.792e8,3.565e8,3.166e8,2.836e8];

channel = struct('label', '', ...
                 'chanpwr', zeros(1,4096));
ch=repmat(channel,1,48);
b = struct('pwrmix', zeros(1,4096), ...
           'divl', zeros(1,4096), ...
           'divu', zeros(1,4096), ...
           'chan', ch);

for i=1:48;
  b.chan(i).chanpwr=kappa(i)*a.chan(i).pwr;
end;

b.pwrmix = 0.0;
b.divl = 0.0;
b.divu = 0.0;

% Calculate power radiated from lower divertor region
for i = 25:31;
  b.divl = b.divl + b.chan(i).chanpwr;
end;

% Calculate power radiated from upper divertor region
for i = 22:24;
  b.divu = b.divu + b.chan(i).chanpwr;
end;

% Calculate total radiated power (based on Tony Leonard's IDL code)
for i = 1:21;
  b.pwrmix=b.pwrmix+b.chan(i).chanpwr;
end;
for i = 6:12;
  b.pwrmix = b.pwrmix - kappa(i)*b.divl/7.0/kappa(i+43);
end;
b.pwrmix = b.pwrmix + b.divu + b.divl;

return;
end
