function [F] = bubfilt(z,sigma)
N = length(z);
V = zeros(N);
V(:,1) = 0.5*(1-z);
V(:,2) = 0.5*(1+z);
for j=3:N
    V(:,j)=LegendreP([zeros(j-3,1);-1;0;1], z);
end

if(nargin==1)
sigma = ones(N,1);
p = 1/16; p=0;
cut = max(1,ceil(N*p));
alpha = (5E-1/2)*max(1/2,1-p*cut);
for k=1:cut
    j = N+1-k;
    sigma(j) = 1-alpha*(((cut+1-k)/cut)^2);
end
end
F = V*diag(sigma)/V;
end