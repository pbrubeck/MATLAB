function []=DirichletSquare(n)
% Solves 2D boundary-value problem specified over the region [-1,1]^2
[D,x]=chebD(n); y=x';
[xx, yy]=ndgrid(x);
D2=D*D;
% Boundary conditions
b1=[0.2*sin(3*pi*y); 0*y];
b2=[(x<0).*sin(pi*x).^4, 0*x];
uu=zeros(n);
uu([1 end],:)=b1;
uu(:,[1 end])=b2;

% Poisson
F=zeros(n);
RHS=F-D2(:,[1 end])*b1-b2*D2(:,[1 end])';
uu(2:end-1, 2:end-1)=sylvester(D2(2:end-1, 2:end-1), D2(2:end-1, 2:end-1)', RHS(2:end-1, 2:end-1));

figure(1);
surf(xx,yy,uu); colormap(jet(256));
camlight; shading interp; axis square;
end