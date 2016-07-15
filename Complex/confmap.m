function [] = confmap(N)
x=linspace(-pi,pi,N);
y=linspace(0,0.8,N);
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=sin(zz);

mesh(real(ww), imag(ww), zeros(N));
axis equal; view(2); colormap([0 0 0]);
end