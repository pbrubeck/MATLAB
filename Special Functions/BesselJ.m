function J=BesselJ(a, x)
t=1;
n=30;
c=ones(1,n);
c(1)=exp(-t)/Gamma(a+1);
for i=1:n-1
    c(i+1)=c(i)*t/(i+a);
end
P=c*Laguerre(n, a);
J=Horner(P, x.^2/(4*t)).*(x/2).^a;
end
