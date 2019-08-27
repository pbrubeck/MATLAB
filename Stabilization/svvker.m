function [Q] = svvker(x)
n=length(x);
V=VandermondeLeg(x);
cut=floor(2*sqrt(n+1)-1);

k=cut+1:n;
sigma=zeros(n,1);
sigma(k) = exp(-((k-n)./(k-cut)).^2);
Q=V'\diag(sigma)/V;
end