function [] = bham2(m,n)
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
r0=L;           % Numerical window
m1=1;           % Mass of BH1
m2=1;           % Mass of BH2
z1= 4;          % Position of BH1
z2=-4;          % Position of BH2
s1= 0.5;        % Spin of BH1
s2= 0.5;        % Spin of BH2

% Differential operators
[Dx,x]=chebD(2*m);
A1=Dx*Dx+diag(1./x)*Dx;
[A1,Dx,x]=radial(A1,Dx,x);
[Dy,y]=chebD(n); y=y';
A2=Dy*Dy;

% Boundary conditions
a=[1,1;1,1];
b=[1,-1;1,-1];
E1=eye(m);
E2=eye(n);
B1=a(1,1)*E1(1,:)+b(1,1)*Dx(1,:);
B2=diag(a(2,:))*E2([1,end],:)+diag(b(2,:))*Dy([1,end],:);
b1=0*y;
b2=[0*x,0*x];

% Coordinate mapping
r=r0*x;
z=r0*y;
[rr,zz]=ndgrid(r,z);
r1=hypot(rr,zz-z1);
r2=hypot(rr,zz-z2);

psi=1+(m1/2)./r1+(m2/2)./r2;
E=r0^2*9/4*(rr.^2).*(s1^2./r1.^8+s2^2./r2.^8+2*s1*s2*(rr.^2+(zz-z1).*(zz-2))./(r1.*r2).^5);
F=zeros(m,n);

[green,ps,kd]=elliptic(A1,A2,B1,B2,1,[1,n]);
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

uu=1+um{1};

figure(1);
surf([rr;-rr(end:-1:1,:)],zz([1:end,end:-1:1],:),uu([1:end,end:-1:1],:));

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