function [] = bipolarMap(R1,R2,z1,z2,N)

L=abs(z1-z2)/2;

N(1:2)=N;
a=sqrt((R1^2-R2^2)^2-8*L^2*(R1^2+R2^2)+16*L^4)/(4*L);
x1=-asinh(a/R2);
x2= asinh(a/R1);

[~,x]=chebD(N(1));
x=x1+(x2-x1)/2*(x+1);
%x=linspace(x1,x2,N(1));

y=linspace(0,2*pi,N(2));
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=z1+(z1-z2)/(2*L)*a*(coth(zz/2)-coth(x2));

mesh(real(ww), imag(ww), 0*xx);
colormap(gray(1)); axis equal;

d=2*(L+R1+R2);
z0=(z1+z2)/2;
xlim(real(z0)+[-d,d]);
ylim(imag(z0)+[-d,d]);

view(2);
end