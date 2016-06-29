function [] = ngon(n, N)
% Solves Poisson's equation on the L-membrane via conformal mapping
p=polygon(exp(2*pi*1i*(0:n-1)/n));
f=diskmap(p);
f=center(f,0);

[Dtt, th] = periodicD2(N);
[Dr,r]=chebD(2*N);
Drr=(diag(r)*Dr)^2;
R2=diag(r.^2);

zz=r*exp(1i*th);
xx=real(zz);
yy=imag(zz);

ww=f(zz);
uu=real(ww);
vv=imag(ww);
df=evaldiff(f,zz);
J=abs(df).^2;

B=-1+zeros(2*N,N);
RHS=J.*(R2*B);

phi=zeros(2*N,N);
phi(2:end-1,:)=sylvester(Drr(2:end-1,2:end-1),Dtt',RHS(2:end-1,:));

phi=[phi(1:N,end), phi(1:N,:)];
uu=[uu(1:N,end), uu(1:N,:)];
vv=[vv(1:N,end), vv(1:N,:)];

figure(1);
surf(uu,vv,phi);
shading interp;
colormap(jet(256));
zrange=max(phi(:))-min(phi(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 4*zrange/hypot(xrange,yrange)]);
view(0,90);
end