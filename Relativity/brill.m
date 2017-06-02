function [] = brill(m,n)
if nargin<2
    n=m;
end

% Simulation parameters
L=8;
r0=sqrt(2)*L;
A0=1;

[Dx,x]=chebD(2*m);
A1=diag(x.^2)*Dx*Dx+diag(2*x)*Dx;
[A1,Dx,x]=radial(A1,Dx,x);

[Dy,y]=chebD(n); y=y';
A2=diag(1-y.^2)*Dy*Dy-diag(2*y)*Dy;

a=[1,1;0,0];
b=[1,1;1,1];

% Imposition of boundary conditions
E1=eye(m);
E2=eye(n);
B1=a(1,1)*E1(1,:)+b(1,1)*Dx(1,:);
B2=diag(a(2,:))*E2([1,end],:)+diag(b(2,:))*Dy([1,end],:);

b1=0*y+1;
b2=[0*x,0*x];

% Coordinate mapping
r=r0*x;
th=acos(y);
rho=r*sin(th);
z=r*cos(th);

qq=A0*(rho.^2).*exp(-(rho/1).^2-(z/1).^2);
C=1/4*((diag(x.^2)*Dx*Dx+diag(x)*Dx)*qq+qq*(diag(1-y.^2)*Dy*Dy-diag(y)*Dy)');
F=zeros(m,n);

uu=elliptic(A1,A2,B1,B2,C,F,b1,b2,1,[1,n],20,1e-10,-1);

figure(1);
surf(kron([-1,1],rho), z(:,[end:-1:1,1:end]), uu(:,[end:-1:1,1:end]));

colormap(jet(256));
colorbar;
shading interp; 
%camlight; 
axis square;
xlim([-L,L]);
ylim([-L,L]);
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$\rho$');
ylabel('$z$');
title('$\psi(\rho,z)$');
view(2);

% [rr,zz]=ndgrid(linspace(-L,L,m));
% xq=hypot(rr,zz)/r0;
% yq=zz./hypot(rr,zz);
% [xx,yy]=ndgrid(x,y);
% psi=interp2(xx',yy',uu',xq',yq')';
% 
% figure(2);
% surf(rr, zz, psi);
% 
% colormap(jet(256));
% colorbar;
% shading interp; 
% %camlight; 
% axis square;
% xlim([-L,L]);
% ylim([-L,L]);
% set(gcf,'DefaultTextInterpreter','latex');
% set(gca,'TickLabelInterpreter','latex','fontsize',14);
% xlabel('$\rho$');
% ylabel('$z$');
% title('$\psi(\rho,z)$');
% view(2);

end