function [ww] = bipolarMap(R1,R2,L,N)
N(1:2)=N;
a=sqrt((R1^2-R2^2)^2-8*L^2*(R1^2+R2^2)+16*L^4)/(4*L);
x1=asinh(a/R1);
x2=asinh(a/R2);

x=linspace(-x1,x2,N(1));
y=linspace(0,2*pi,N(2));
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=a*coth(zz/2);%-a*coth(x2)+L;
end