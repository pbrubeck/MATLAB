function [] = DirichletEllipse(N,M)
% Solves Dirichlet BVP over an elliptical region
[Gu,Gv,Duu,Dvv,u,v]=chebLapEll(N,M);
d1=Duu(:,1);
dn=Duu(:,N);
f=1; % focal distance/2
xx=f*cosh(u)*cos(v);
yy=f*sinh(u)*sin(v);

% Boundary conditions
g=1*cos(5*v);
F=zeros(N,M);
RHS=f^2*(Gu*F+F*Gv)-(d1+dn)*g;

Psi=zeros(N,M);
Psi(1,:)=g;
Psi(2:end-1,:)=sylvester(Duu(2:end-1,2:end-1), Dvv', RHS(2:end-1,:));

half=floor((N+1)/2);
xx=xx(1:half,:);
yy=yy(1:half,:);
Psi=Psi(1:half,:);

xx=[xx(:,end), xx];
yy=[yy(:,end), yy];
zz=real([Psi(:,end), Psi]);
figure(1); surfl(xx,yy,zz,'light'); 
shading interp; colormap(jet(256));
zrange=max(zz(:))-min(zz(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 zrange/hypot(xrange,yrange)]);
end