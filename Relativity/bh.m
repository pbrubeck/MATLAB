function [] = bh(m,n)
% Black-hole, excision method
if nargin<2
    n=m;
end

% Homotopy Analysis Method
its=30;
tol=1e-6;
h=-1;

% Simulation parameters
L=20;   % Plot window
a0=1;   % Throat
p=0.5;  % Momentum

% Differential operators
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
A1=diag((x+1).^2)*Dx*Dx;
A2=diag(1-y.^2)*Dy*Dy-diag(2*y)*Dy;
E=3/256*(p/a0)^2*((x+1).^2.*(3-2*x-x.^2).^2);
psi=1;
F=zeros(m,n);

% Coordinate mapping
r=(2*a0)./(x+1);
th=acos(y);
rho=r*sin(th);
z=r*cos(th);

% Boundary conditions
a=[-1/4,1;0,0];
b=[1,0;1,1];
b1=[0*y; 0*y];
b2=[0*x, 0*x];

E1=eye(m);
E2=eye(n);
B1=diag(a(1,:))*E1([1,end],:)+diag(b(1,:))*Dx([1,end],:);
B2=diag(a(2,:))*E2([1,end],:)+diag(b(2,:))*Dy([1,end],:);

% Solution
[green,ps,kd]=elliptic(A1,A2,B1,B2,[1,m],[1,n]);
eqn=@(uu,F) kd(A1*uu+uu*A2'+E.*(psi+uu).^(-7)-F);

% HAM nonlinear functions
R1=@(um) (psi+um{1}).^(-7);
R2=@(um) -7*(psi+um{1}).^(-8).*(um{2});
R3=@(um) -7*(psi+um{1}).^(-9).*((psi+um{1}).*um{3}-4*um{2}.^2);
R4=@(um) -7*(psi+um{1}).^(-10).*((psi+um{1}).^2.*um{4}-8*(psi+um{1}).*um{2}.*um{3}+12*um{2}.^3);

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
uu=psi+um{1};

% Plot
[~,mm]=max(r>=sqrt(2)*L);
figure(1);
surf(kron([-1,1],rho(1:mm,:)),z(1:mm,[end:-1:1,1:end]), uu(1:mm,[end:-1:1,1:end]));
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