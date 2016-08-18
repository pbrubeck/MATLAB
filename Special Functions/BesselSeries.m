function [y] = BesselSeries(a, x)
% Bessel Series: sum(a(n)*J(n,x))
n=length(a)-1;
y=zeros(size(x));
yy=zeros(size(x));
for k=0:n-1
    temp=y;
    y=-(yy-2*k.*y./x-a(k+1));
    yy=temp;
end
jn=besselj(n,x);
jn1=besselj(n-1,x);
y=a(n+1)*jn+jn1.*y+jn.*yy;
end