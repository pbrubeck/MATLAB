function [] = flt(n)
% Fast Legendre Transform, based on Alpert-Rokhlin

% Initialization
k=32;
s=4*k;

% Step 1
r=0:k-1;
t=(1-cos((r+0.5)*pi/k))/2;
tt=[t/2, (1+t)/2];

% Step 2
T=repmat(t,[k,1]);
den=prod(T-T'+eye(k));
l=0:s-1;
U=repmat(l/s,[k,1])-repmat(t(:),[1,s]);
U=bsxfun(@rdivide, bsxfun(@rdivide, prod(U), U), den(:));

%Step 3
h=log2(n/s)-1;

end

