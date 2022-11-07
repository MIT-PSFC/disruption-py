function y = hp(x, t, fLow)
%Highpass function - implement a lowpass filter through Fourier transform.
% Gifted to Alex Tinguely by Ted Golfinopolous 04-29-16
%USAGE:
%
% y = hp( x, t, fLow )
%
%  x = input signal
%  t = time base of signal OR sampling time of signal
%  fLow = lower bound of frequency pass band [Hz]
%
% Note: fLow must be less than the Nyquist frequency.
%
% y = filtered signal
%
%Ted Golfinopoulos, 13 April 2010
%See also myFft

if(length(t)>1)
    ts=t(2)-t(1); %Sampling time.
else
    ts=t;
end

%Error if frequency band extends past Nyquist frequency.
if( fLow > 1/(2*ts) )
    error('work/hp.m:bandExceedsNyquist', 'Lower frequency bound must be below the Nyquist frequency.');
end

%Calculate one-sided Fourier transform (two-sided transform, in my
%implementation,  argument causes myFft to return two-sided transform.
[X, f] = myFft(x, t);

%Implement high-pass filter.
%Frequencies outside of pass band are nullified.
X( f<fLow | (f > f(end)-fLow)  )=0;

%Transform back to time domain.
y = myIfft(X,f);