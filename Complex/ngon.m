function [] = ngon(n, N)
% Solves Poisson's equation on the L-membrane via conformal mapping
p=polygon(exp(2*pi*1i*(0:n-1)/n));
f=diskmap(p);
f=center(f,0);

[R2,Drr,Dtt,r,th]=chebLapPol(N,N);
zz=r*exp(1i*th);

ww=f(zz);
uu=real(ww);
vv=imag(ww);
df=evaldiff(f,zz);
J=abs(df).^2;

B=-1+zeros(N,N);
RHS=(R2*J).*(B);

phi=zeros(N,N);
phi(2:end,:)=sylvester(Drr(2:end,2:end),Dtt',RHS(2:end,:));

phi=[phi(1:N,end), phi(1:N,:)];
uu=[uu(1:N,end), uu(1:N,:)];
vv=[vv(1:N,end), vv(1:N,:)];

figure(1);
mesh(uu,vv,phi);
shading interp;
colormap(jet(256));
zrange=max(phi(:))-min(phi(:));
xrange=max(uu(:))-min(uu(:));
yrange=max(vv(:))-min(vv(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
view(0,90);
end