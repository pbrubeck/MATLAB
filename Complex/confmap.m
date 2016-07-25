function [] = confmap(N)
x=linspace(0,0.5,N);
y=linspace(0,2*pi,N);
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=cosh(zz);

mesh(real(ww), imag(ww), zeros(N));
axis equal; view(2); colormap([0 0 0]);
end