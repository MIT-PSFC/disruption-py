function [y, t] = myIfft(Y, f, varargin)
%Given discrete Fourier spectrum, Y, for frequencies, f, returns
%time-domain signal, y, with time vector, t.  Inverts the fft function,
%accounting for normalization by time span of that function.
% Gifted to Alex Tinguely by Ted Golfinopolous 04-29-16
%USAGE:
%
% [y, t] = myIfft(Y, f[, twoSidedFlag])
% [y, t] = myIfft(Y, df[, twoSidedFlag])
%
%  t(1) = 0
%
% If the second argument contains only one element (scalar), it is taken to
% be the difference between discrete frequency points, or one over the time
% span.
%
% If twoSidedFlag input is included and is set to true, spectrum must be
% input as two-sided from -fN to fN.  Note that the current implementation
% of myFft destroys some information in the spectrum when it produces a
% two-sided spectrum, so that 
%
%Ted Golfinopoulos, Thursday, 11 February 2010
%See also myFft, fft, ifft

%Reverse two-sided transform if needed.
if(nargin>2)
    if(varargin{1})
        mid = ceil(length(Y)/2);
        Y(mid)=Y(mid)*0.5;
        Y = [Y(mid:end); Y(1:mid)];
    end
end

if(length(f)==1)
    DeltaF=f;% argument, f, is df.
else
    DeltaF=f(2)-f(1);
end

DeltaT=1/((length(Y)-1)*DeltaF);

y = ifft( Y / DeltaT );

%Only calculate if requested.
if(nargout>1)
    t = linspace(0, DeltaT, length(y));
end