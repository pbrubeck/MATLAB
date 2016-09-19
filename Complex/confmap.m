function [] = confmap(a,b,N)
f=sqrt(a^2-b^2);
xi0=acosh(a/f);

x=linspace(0,xi0,N);
y=linspace(0,2*pi,N);
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=f*cosh(zz);

mesh(real(ww), imag(ww), zeros(N));
axis equal; view(2); colormap([0 0 0]);
end