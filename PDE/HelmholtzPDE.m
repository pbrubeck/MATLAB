function []=HelmholtzPDE(n, k)
% Solves the Helmholtz equation in two dimensions using Chebyshev methods
[D, x]=chebD(n); y=x';
[xx, yy]=ndgrid(x);
D2=D^2+(k^2/2)*eye(n); 
d1=D2(:,1);
dn=D2(:,end);
% Boundary conditions
ua=0*y;
ub=0.2*sin(3*pi*y);
uc=0*x;
ud=(x<0).*sin(pi*x).^4;
whos

F=0*exp(-10*((xx-.5).^2+(yy-1).^2));
RHS=F-dn*ua-d1*ub-uc*dn'-ud*d1';
uu=zeros(n);
uu(n,:)=ua; uu(1,:)=ub;
uu(:,n)=uc; uu(:,1)=ud;
uu(2:end-1, 2:end-1)=sylvester(D2(2:end-1, 2:end-1), D2(2:end-1, 2:end-1)', RHS(2:end-1, 2:end-1));

figure(1);
surf(xx,yy,uu); colormap(jet(256));
camlight; shading interp; axis square;
end