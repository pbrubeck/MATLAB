function [] = bham(m,n)
% Black-hole with angular momentum, puncture method
if nargin<2
    n=m;
end

% Homotopy Analysis Method
maxits=30;
tol=3e-10;
h=-1;

% Simulation parameters
L=20;           % Plot window
r0=sqrt(2)*L;   % Numerical window
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
E=9/4*J^2*(r.^-4)*(1-y.^2);
F=zeros(m,n);

[green,ps,kd]=elliptic(A1,A2,B1,B2,1,[1,n]);
eqn=@(uu,F) kd(A1*uu+uu*A2'+E.*(psi+uu).^(-7)-F);

% HAM nonlinear functions
R{1}=@(um) (psi+um{1}).^(-7);
R{2}=@(um) -7*(psi+um{1}).^(-8).*(um{2});
R{3}=@(um) -7*(psi+um{1}).^(-9).*((psi+um{1}).*um{3}-4*um{2}.^2);
R{4}=@(um) -7*(psi+um{1}).^(-10).*((psi+um{1}).^2.*um{4}-8*(psi+um{1}).*um{2}.*um{3}+12*um{2}.^3);
R{5}=@(um) -7*(psi+um{1}).^(-11).*(-30*um{2}.^4+36*(psi+um{1}).*um{2}.^2.*um{3}+...
    -8*(psi+um{1}).^2.*um{2}.*um{4}+(psi+um{1}).^2.*(-4*um{3}.^2+(psi+um{1}).*um{5}));

ub=ps(b1,b2);
uu=ub-green(eqn(ub,F));

its=0;
err=1;
um=cell(length(R)+1,1);
while err>tol && its<maxits
    um{1}=uu;
    for k=1:length(R)
        um{k+1}=(k>1)*um{k}+h*green(kd(A1*um{k}+um{k}*A2'+E.*R{k}(um)));
    end
    uu=sum(cat(3,um{:}),3);
    err=norm(1-um{1}./uu,'inf');
    its=its+1;
end

display(its);
display(err);
uu=uu+1;

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