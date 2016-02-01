function [y] = BesselS(a, x)
% Bessel Series: sum(a(n)*J(n,x))
n=length(a);
y=(n>1)*a(n); 
yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    y=a(k+1)+2*k*y./x-yy;
    yy=temp;
end
j0=besselj(0,x);
j1=besselj(1,x);
y=a(1)*j0+y.*j1-yy.*j0;
end