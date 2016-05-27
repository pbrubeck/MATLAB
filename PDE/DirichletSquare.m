function []=DirichletSquare(n)
% Solves 2D boundary-value problem specified over the region [-1,1]^2
[D,x]=chebD(n); y=x';
[yy, xx]=meshgrid(x);
D2=D*D; d1=D2(:,1); dn=D2(:,end);
% Boundary conditions
ua=0*y;
ub=0.2*sin(3*pi*y);
uc=0*x;
ud=(x<0).*sin(pi*x).^4;

% Poisson
F=zeros(n);
rhs=F-dn*ua-d1*ub-uc*dn'-ud*d1';

tic;
uu=zeros(n);
uu(n,:)=ua; uu(1,:)=ub;
uu(:,n)=uc; uu(:,1)=ud;
uu(2:end-1, 2:end-1)=lyap(D2(2:end-1, 2:end-1), -rhs(2:end-1, 2:end-1));
toc
disp(norm(D2*uu+uu*D2'-F, 'fro'));

figure(1); surfl(xx,yy,uu,'light'); 
shading interp; colormap(jet(256));
end