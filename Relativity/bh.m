function [] = bh(m,n)
% Black-hole
if nargin<2
    n=m;
end

% Simulation parameters
L=8;  % Window
a0=1; % Throat
p=2;  % Momentum

% Differential operators
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
A1=diag((x+1).^2)*Dx*Dx;
A2=diag(1-y.^2)*Dy*Dy-diag(2*y)*Dy;
R=3/256*(p/a0)^2*((x+1).^2.*(3-2*x-x.^2).^2);
C=@(uu) R.*(uu.^(-7));
F=zeros(m,n);

% Coordinate mapping
r=(2*a0)./(x+1);
th=acos(y);
rho=r*sin(th);
z=r*cos(th);

% Boundary conditions
a=[-1/4,1;0,0];
b=[1,0;1,1];
b1=[0*y; 0*y+1];
b2=[0*x, 0*x];

E1=eye(m);
E2=eye(n);
B1=diag(a(1,:))*E1([1,end],:)+diag(b(1,:))*Dx([1,end],:);
B2=diag(a(2,:))*E2([1,end],:)+diag(b(2,:))*Dy([1,end],:);

% Solution
[uu]=elliptic(A1,A2,B1,B2,C,F,b1,b2,[1,m],[1,n]);

% Plot
[~,mm]=max(r>=sqrt(2)*L);
figure(1);
surf(kron([-1,1],rho(1:mm,:)),z(1:mm,[end:-1:1,1:end]), uu(1:mm,[end:-1:1,1:end]));
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
end