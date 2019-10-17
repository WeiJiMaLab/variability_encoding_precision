function Y = circ_vmrnd(theta, kappa, n)
% alpha = circ_vmrnd(theta, kappa, n)
%   Simulates n random angles from a von Mises distribution, with preferred 
%   direction thetahat and concentration parameter kappa.
%
%   Input:
%     theta   - preferred direction
%     kappa   - width
%     n       - number of samples (default: 1)
%
%   Output:
%     alpha   - samples from von Mises distribution
%
%   Theta and kappa can be scalars, vectors, or matrices. If both or not
%   scalar, their dimensions should be identical and n should be 1. 
%
%   Examples
%    Y = circ_vmrnd(0,5)            - draw a sample from VM(0,5)
%    Y = circ_vmrnd(0,5,100)        - draw 100 samples from VM(0,5)
%    Y = circ_vmrnd(0,[5 10],100)   - draw 100 samples from VM(0,5) and 100 from VM(0,10)
%    Y = circ_vmrnd([-pi/2 pi/2],5) - draw one sample from VM(-pi/2,5) and one from VM(pi/2,5)
%
%   References:
%     Statistical analysis of circular data, Fisher, sec. 3.3.6, p. 49
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens and Marc J. Velasco, 2009
% Modified by Ronald van den Berg, 2011 

% input checking
if ~exist('n','var') 
    n=1;
end
if n==0
    Y = [];
    return
end

% handle all cases
if numel(kappa)==1 && numel(theta)==1     % both inputs are scalar
    input_dims = size(kappa);
    kappa=repmat(kappa,1,n);
    theta=repmat(theta,1,n);
elseif numel(kappa)==1 && numel(theta)>1  % kappa is scalar, theta is matrix
    input_dims = size(theta);
    theta = theta(:)';
    theta = repmat(theta,1,n);    
    kappa = ones(size(theta))*kappa;
elseif numel(kappa)>1 && numel(theta)==1  % kappa is matrix, theta is scalar
    input_dims = size(kappa);
    kappa = kappa(:)';
    kappa = repmat(kappa,1,n);    
    theta = ones(size(kappa))*theta;
elseif numel(kappa)>1 && numel(theta)>1   % both inputs are matrices
    if n>1
        error('Can only have n>1 when theta and kappa are scalars or vectors');        
    end
    if ~isequal(size(theta),size(kappa))
        error('Invalid input dimensions. If both theta and kappa is a vector or matrix, their dimensions should be the same');        
    end
    input_dims = size(kappa);
    kappa = kappa(:)';
    theta = theta(:)';
end

% use code from original circ_vmrnd to draw samples
a = 1 + sqrt((1+4*kappa.^2));
b = (a - sqrt(2*a))./(2*kappa);
r = (1 + b.^2)./(2*b);
valid = zeros(1,length(kappa));
z = zeros(size(kappa));
f = zeros(size(kappa));
c = zeros(size(kappa));
while ~all(valid)
    u(:,~valid) = rand(3,sum(~valid));    
    z(~valid) = cos(pi*u(1,~valid));
    f(~valid) = (1+r(~valid).*z(~valid))./(r(~valid)+z(~valid));
    c(~valid) = kappa(~valid).*(r(~valid)-f(~valid));       
    valid = u(2,:) < c .* (2-c) | ~(log(c)-log(u(2,:)) + 1 -c < 0);               
end
Y = theta + sign(u(3,:) - 0.5) .* acos(f);
Y = angle(exp(1i*Y));

% if kappa is very small, draw from a uniform distribution (the above method gives NaN for very small kappa's)
Y(kappa<1e-6) = 2*pi*rand(sum(kappa<1e-6),1)-pi;

% reshape back to original dimensions of input (add a dimension when n>1)
if n>1
    Y = reshape(Y,[max(input_dims) n]);
else
    Y = reshape(Y,input_dims);
end
