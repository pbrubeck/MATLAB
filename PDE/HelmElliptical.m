function [] = HelmElliptical(N,m)
f=1;
[Gu,Gv,Duu,Dvv,u,v]=chebLapEll(2*N-1,N);
Duu=Duu(2:end-1,2:end-1);
Gu=Gu(2:end-1,2:end-1);

a=m^2;
lam=-2*a;
Vu=cos(pi*u(2:end-1));
Vv=cos(m*v(:));
[Vu,Vv,lam,a]=mpep(Duu,Dvv,Gu,Gv,Vu,Vv,lam,a);
q=-lam/4;

rmf=zeros(N,1);
rmf(2:end,:)=bsxfun(@times, conj(Vu(N-1,:))./abs(Vu(N-1,:)), Vu(1:N-1,:));
amf=bsxfun(@times, conj(Vv(1,:))./abs(Vv(1,:)), Vv);
u=u(1:N);

display(a);
display(q);
figure(1); plot(u,rmf);
figure(2); plot(v,amf);

xx=f*cosh(u)*cos([0,v]);
yy=f*sinh(u)*sin([0,v]);
zz=rmf*amf.';
zz=[zz(:,end), zz];
figure(3);
surf(xx,yy,zz); 
shading interp;
axis off;
colormap(jet(256));
zrange=max(zz(:))-min(zz(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 zrange/hypot(xrange,yrange)]);
view(0,90);
end