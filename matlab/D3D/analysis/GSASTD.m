function Processed = GSASTD(x,y,DerivativeMode,width,type,ends,SlewRate)

% 'DerivativeMode' determines the derivative order:
% -------------------------------------------------
%       0 -> smoothing
%       1 -> first derivative

% 'width' determines the smooth width:
% ------------------------------------

% 'type' determines the smooth mode:
% ----------------------------------
%       0 -> no smoothing.
%       1 -> rectangular (sliding-average or boxcar)
%       2 -> triangular (2 passes of sliding-average)
%       3 -> pseudo-Gaussian (3 passes of sliding-average)

% 'ends' controls how the "ends" of the signal are handled.
%       0 -> ends are "zeroed"
%       1 -> the ends are smoothed with progressively smaller smooths the closer to the end.

if SlewRate
    for n=1:length(y)-1
        if y(n+1)-y(n)>SlewRate,y(n+1)=y(n)+SlewRate;end
        if y(n+1)-y(n)<-SlewRate,y(n+1)=y(n)-SlewRate;end
    end
end % SlewRate

% ----------------------------------------------------------------------
if type==0,width=1;end
switch DerivativeMode
    case 0
        Processed=fastsmooth(y,width,type,ends);
    case 1
        Processed=fastsmooth(deriv(x,y),width,type,ends);
end

Processed=reshape(Processed,size(x));

% ----------------------------------------------------------------------
function d=deriv(x,y)
n=length(y);
d=zeros(size(y));
d(1)=(y(2)-y(1))./(x(2)-x(1));
d(n)=(y(n)-y(n-1))./(x(n)-x(n-1));
for j = 2:n-1
    d(j)=(y(j+1)-y(j-1)) ./ (1.*(x(j+1)-x(j-1)));
end

% ----------------------------------------------------------------------
function SmoothY=fastsmooth(Y,w,type,ends)
if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
switch type
    case 0
        SmoothY=sa(Y,w,ends);
    case 1
        SmoothY=sa(Y,w,ends);
    case 2
        SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
        SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
end
function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w
    s(k+halfw-1)=SumPoints;
    SumPoints=SumPoints-Y(k);
    SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;

% ----------------------------------------------------------------------
if ends==1
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint
        SmoothY(k)=mean(Y(1:(2*k-1)));
        SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
end% ----------------------------------------------------------------------
