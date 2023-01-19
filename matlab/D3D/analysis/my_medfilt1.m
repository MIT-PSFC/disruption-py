function y = my_medfilt1(x,n,method,blksz,DIM) 

if nargin < 2, n = []; end
if nargin < 3, method = 'nearest'; end
if nargin < 4, blksz = []; end
if nargin < 5, DIM = []; end

% Ends treatments
x = x(:);
x = [repmat(mean(x(1:2)),round(n/2),1)+randn(round(n/2),1);x;repmat(mean(x(end-1:end)),round(n/2),1)+randn(round(n/2),1)];

% Check if the input arguments are valid
if isempty(n)
  n = 3;
end

if ~isempty(DIM) && DIM > ndims(x)
	error(message('signal:medfilt1:InvalidDimensions'))
end

% Reshape x into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first (row) dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

% Verify that the block size is valid.
siz = size(x);
if isempty(blksz)
	blksz = siz(1); % siz(1) is the number of rows of x (default)
else
	blksz = blksz(:);
end

% Initialize y with the correct dimension
y = zeros(siz); 

% Call medfilt1D (vector)
for i = 1:prod(siz(2:end))
	y(:,i) = medfilt1D(x(:,i),n,blksz);
end

% Convert y to the original shape of x
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end

% Ends treatment
% y = interp1(2:numel(y),y(2:end),1:numel(y),method,'extrap');
y = y(:);
y_check = [1;diff(y)];
y_base = (1:numel(y_check))';
y = interp1(y_base(y_check~=0),y(y_check~=0),y_base,method,'extrap');
y = y(round(n/2)+1:end-round(n/2));

%-------------------------------------------------------------------
%                       Local Function
%-------------------------------------------------------------------
function y = medfilt1D(x,n,blksz)
%MEDFILT1D  One dimensional median filter.
%
% Inputs:
%   x     - vector
%   n     - order of the filter
%   blksz - block size

nx = length(x);
if rem(n,2)~=1    % n even
    m = n/2;
else
    m = (n-1)/2;
end
X = [zeros(m,1); x; zeros(m,1)];
y = zeros(nx,1);

% Work in chunks to save memory
indr = (0:n-1)';
indc = 1:nx;
for i=1:blksz:nx
    ind = indc(ones(1,n),i:min(i+blksz-1,nx)) + ...
          indr(:,ones(1,min(i+blksz-1,nx)-i+1));
    xx = reshape(X(ind),n,min(i+blksz-1,nx)-i+1);
    y(i:min(i+blksz-1,nx)) = median(xx,1);
end

