function [Y, f] = myFft(y, t, varargin)
%Returns fast fourier transform of y(t), as well as frequency space.
% Gifted to Alex Tinguely by Ted Golfinopolous 04-29-16
%Usage:
%
% [Y, f] = myFft(y, t[,twoSided])
% [Y, f] = myFft(y, ts[,twoSided])
%
% y and t are vectors of the same length, and t should have uniform
% spacing and be always increasing (or decreasing).
%
% twoSided is optional argument.  If true, make Y on a two-sided scale,
% with f from -fNyquist to fNyquist.  Otherwise, f from 0 to 2*fNyquist
% (default).
%
% If the second argument is a scalar (single element), it is taken to be
% the time interval between samples, or one over the sampling frequency.
%
%
%Ted Golfinopoulos, 21 September 2009 (Fall Equinox)
%Last updated Thursday: 11 February 2010: note that will not work with ifft
%because output, Y, is normalized by max(t)-min(t).  Use myIfft instead.
%See also myIfft, fft, ifft

if(length(t)==1)
    df = 1/t; %Sampling frequency given as argument instead of time vector.
    deltaT = length(y)*t; % Total time elapsed.
else
    dt = t(2) - t(1);
    df = 1.E0 / dt;
    deltaT = max(t) - min(t);
end

% Y = fft(y) / deltaT; %Changed 17 July 2012
Y = fft(y) / df; %Multiply by sampling time to get Fourier transform amplitude rather than fft amplitude that's not normalized.

%Make Y a row-vector.
if(size(Y,2) > size(Y,1) ),
    Y = Y.';
end

f = df * linspace(0,1, length(Y)).'; %Make a row-vector.

%Give two-sided transform.
if(nargin>2)
    if(varargin{1})
        mid = ceil(length(Y)/2);
%         Y = [Y(mid+1:end-1); Y(end)+Y(1); Y(2:mid)];
%         f = [f(mid+1:end)-f(end); f(2:mid)];
        Y=[Y(mid+2:end); Y(1:mid+1)];
        f=[f(mid+2:end)-f(end)-(f(2)-f(1)); f(1:mid+1)];
    end
end

