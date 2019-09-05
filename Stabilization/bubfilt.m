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
p = 0.0;
cut = max(1,ceil(N*p));
alpha = 0.01;
for k=1:cut
    w = alpha*(((cut+1-k)/cut)^2);
    j = N+1-k;
    sigma(j) = 1-w;
end
end

F = V*diag(sigma)/V;
end