function [] = brillbh(m,n)
% Brill wave plus black hole
if nargin<2
    n=m;
end

% Simulation parameters
L=8;    % Window
a0=1;   % Throat
A0=1;   % Amplitude
s0=1; 
eta0=1;
n0=2;

% Differential operators
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
A1=diag((x+1).^2)*Dx*Dx;
A2=diag(1-y.^2)*Dy*Dy-diag(2*y)*Dy;

% Coordinate mapping
r=(2*a0)./(x+1);
th=acos(y);
rho=r*sin(th);
z=r*cos(th);

eta=log(2./(x+1));
qq=A0*(exp(-((eta+eta0)/s0).^2)+exp(-((eta-eta0)/s0).^2))*((1-y.^2).^(n0/2));
qq(isnan(qq))=0;
C=1/4*((diag((x+1).^2)*Dx*Dx+diag(x+1)*Dx)*qq+qq*(diag(1-y.^2)*Dy*Dy-diag(y)*Dy)');
F=zeros(m,n);

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
[green,ps,kd,sc,gb]=elliptic(A1,A2,B1,B2,[1,m],[1,n]);
afun=@(uu) sc(uu)+kd(C).*uu;
pfun=@(uu) kd(green(uu));

ub=ps(b1,b2);
rhs=kd(F-A1*ub-ub*A2'-C.*ub);
u0=kd(ub+green(rhs));

[uu,res,its]=precond(afun,pfun,rhs,u0,20,1e-10);
uu=gb(uu)+ub;
display(res);
display(its);

% Plot
[~,mm]=max(r>=sqrt(2)*L);
figure(1);
surf(kron([-1,1],rho(1:mm,:)), z(1:mm,[end:-1:1,1:end]), uu(1:mm,[end:-1:1,1:end]));

colormap(jet(256));
colorbar;
shading interp;
camlight; 
axis square;
xlim([-L,L]);
ylim([-L,L]);
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$\rho$');
ylabel('$z$');
title('$\psi(\rho,z)$');
end