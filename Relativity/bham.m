function [] = bham(m,n)
% Black-hole with angular momentum, puncture method
if nargin<2
    n=m;
end

% Homotopy Analysis Method
its=30;
tol=1e-7;
h=-1;

% Simulation parameters
L=20;           % Plot window
r0=sqrt(2)*L;   % Radial window
m1=1;           % Mass
J=0.5;          % Angular momentum

% Differential operators
[Dx,x]=chebD(2*m);
A1=diag(x.^2)*Dx*Dx+diag(2*x)*Dx;
[A1,Dx,x]=radial(A1,Dx,x);
[Dy,y]=chebD(n); y=y';
A2=diag(1-y.^2)*Dy*Dy-diag(2*y)*Dy;

% Boundary conditions
a=[1,1;0,0];
b=[1,1;1,1];
E1=eye(m);
E2=eye(n);
B1=a(1,1)*E1(1,:)+b(1,1)*Dx(1,:);
B2=diag(a(2,:))*E2([1,end],:)+diag(b(2,:))*Dy([1,end],:);
b1=0*y;
b2=[0*x,0*x];

% Coordinate mapping
r=r0*x;
th=acos(y);
rho=r*sin(th);
z=r*cos(th);
rr=rho+1i*z;

psi=1+(m1/2)./abs(rr);
E=9/4*J^2*(r.^-4).*(1-y.^2);
F=zeros(m,n);

[green,ps,kd]=elliptic(A1,A2,B1,B2,1,[1,n]);
eqn=@(uu,F) kd(A1*uu+uu*A2'+E.*(psi+uu).^(-7)-F);

% HAM nonlinear functions
R1=@(um) (psi+um{1}).^(-7);
R2=@(um) -7*(psi+um{1}).^(-8).*(um{2});
R3=@(um) -7/2*(psi+um{1}).^(-9).*(2*(psi+um{1}).*um{3}-8*um{2}.^2);
R4=@(um) -7/6*(psi+um{1}).^(-10).*(6*(psi+um{1}).^2.*um{4}-48*(psi+um{1}).*um{2}.*um{3}+72*um{2}.^3);

ub=ps(b1,b2);
um=cell(4,1);
um{1}=ub-green(eqn(ub,F));

i=0;
res=norm(eqn(um{1},F),'inf');
while res>tol && i<its
    um{2}=h*green(kd(A1*um{1}+um{1}*A2'+E.*R1(um)));
    um{3}=um{2}+h*green(kd(A1*um{2}+um{2}*A2'+E.*R2(um)));
    um{4}=um{3}+h*green(kd(A1*um{3}+um{3}*A2'+E.*R3(um)));
    um{5}=um{4}+h*green(kd(A1*um{4}+um{4}*A2'+E.*R4(um)));
    um{1}=um{1}+um{2}+um{3}+um{4}+um{5};
    res=norm(eqn(um{1},F),'inf');
    i=i+1;
end
display(res);
display(i);

uu=1+um{1};

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
title('$u(\rho,z)$');
view(2);
end