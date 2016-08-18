function [y] = sumj2k(a, s, x)
% Bessel even/odd series a'*besselj(2*k+s, x)
n=length(a);
x2=x.*x;
y=(n>1)*a(n); 
yy=zeros(size(x));
for k=n-2:-1:1
    j=2*k+s;
    temp=y;
    y=a(k+1)+(4*j*(j+1)./x2-2*j/(j-1)).*y-(j+3)/(j+1)*yy;
    yy=temp;
end
js=besselj(s,x);
js2=besselj(s+2,x);
y=a(1)*js+js2.*y-(s+3)/(s+1)*js.*yy;
end
