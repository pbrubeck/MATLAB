function [ww] = ellipticMap(a,b,N)
N(1:2)=N;
f=sqrt(a^2-b^2);
xi0=acosh(a/f);
x=linspace(0,xi0,N(1));
y=linspace(0,2*pi,N(2));
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=f*cosh(zz);
end