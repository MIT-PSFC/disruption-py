function [sigma,mu,A] = gaussian_fit(xdata,ydata)

%---------------------------------------------------------------------
%
% This function outputs the standard deviation, peak center, and peak value
% of a Gaussian fit to the inputs xdata, ydata using a least squares linear
% regression approach. Assuming that y fits to the form
%
%       ydata = 1/(sigma*sqrt(2*pi)) * exp( -1/2*([xdata-mu]/sigma)^2 )
%
% --> ln(ydata) = ln(1/sigma*sqrt(2*pi)) - 1/2*([xdata-mu]/sigma)^2
%
%           = c2 * xdata^2 + c1 * xdata + c0
%
% we can solve the equation Ax = b where 
%
% A = [xdata.^2 xdata ones(size(xdata))] (an nx3 matrix)
% 
% b = ln(ydata) (an nx1 matrix)
%
% c = [c2 c1 c0] (a 3x1 matrix with the coefficients. 
% 
% Reverse engineering the Gaussian parameters from the coefficients of c
% then gives the output.

% Check first to see that input data are nx1 matrices; if not, then
% transpose them to columns

if (size(xdata,2) > 1)
  xdata = xdata';
end;
if (size(ydata,2) > 1)
  ydata = ydata';
end;

% Initialize to null values

sigma = NaN; mu = NaN; A = NaN;

% Check input data for NaN values and remove any corresponding data

bad_indx = union(find(isnan(xdata)),find(isnan(ydata)));
indices = setdiff(1:length(xdata),bad_indx);
xdata = xdata(indices);
ydata = ydata(indices);

% Find least squares solution, gives coeffs x = [c2, c1, c0] 
A = [xdata.^2, xdata, ones(size(xdata))];
b = log(ydata);
x = A\b;

% Output standard deviation (sigma), peak center (mu), and peak value (A) 
mu = -x(2)/x(1)/2;
sigma = sqrt( -1/2/x(1) );
A = exp( x(3)-x(2)^2/4/x(1) );

% Check to see if results make sense. If standard deviation is negative,
% assume there is no good Gaussian fit and output NaN values for each 
% parameter.

if imag(sigma)~=0
    sigma = NaN; mu = NaN; A = NaN;
end
