function [ V ] = VandermondeGeg(a,x)
V=zeros(length(x));
L=pi*2^(1-2*a)*gamma(2*a)/(a*gamma(a)^2);
for j=1:length(x)
    V(:,j)=(1/sqrt(L))*GegenbauerC([zeros(j-1,1);1], a, x);
    L=L*((j-1+a)*(j-1+2*a))/(j*(j+a));
end
end