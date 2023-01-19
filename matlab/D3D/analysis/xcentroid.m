function centroidval = xcentroid(x, y);

centroidval = NaN;
total = 0;

for i = 1:length(y);
  total = total + x(i) * y(i);
end;

centroidval = total / sum(y);
end
