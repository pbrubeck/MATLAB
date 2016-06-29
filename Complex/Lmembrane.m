function [] = Lmembrane( N )
% Solves Poisson's equation on the L-membrane via conformal mapping
p=polygon([1i -1+1i -1-1i 1-1i 1 0]);
f=rectmap(p, [5 1 2 4]);
params=parameters(f);
xmin=min(real(params.prevertex));
xmax=max(real(params.prevertex));
ymin=min(imag(params.prevertex));
ymax=max(imag(params.prevertex));
dx=xmax-xmin;
dy=ymax-ymin;

[D,t]=chebD(N);
[xx,yy]=ndgrid(xmin+dx*(t+1)/2, ymin+dy*(t+1)/2);
zz=xx+1i*yy;

ww=f(zz);
uu=real(ww);
vv=imag(ww);
df=evaldiff(f,zz);
J=abs(df).^2;

D2=D^2;
A=D2(2:end-1,2:end-1);
B=-1+zeros(N,N);
RHS=J.*B;

phi=zeros(N,N);
phi(2:end-1,2:end-1)=sylvester((2/dx)^2*A,(2/dy)^2*A',RHS(2:end-1,2:end-1));

figure(1);
surf(uu,vv,phi);
shading interp;
colormap(jet(256));
zrange=max(phi(:))-min(phi(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 4*zrange/hypot(xrange,yrange)]);
xlim([-1 1]);
ylim([-1 1]);
view(0,90);
end