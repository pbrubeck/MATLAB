function [] = PoissonGalerkin( m, n )
% Solves the Poisson equation using Legedre collocation - weak Galerkin
% spectral method.
if(nargin==1)
    n=m;
end

kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];

% Boundary conditions
a=[1,1;1,1];
b=[0,1;0,1];

[K1,M1,E1,V1,L1,x]=GalerkinGLL(m,a(1,:),b(1,:));
[K2,M2,E2,V2,L2,y]=GalerkinGLL(n,a(2,:),b(2,:));
[Lx,Ly]=ndgrid(L1,L2); LL=Lx+Ly;

function uu=greenF(F,b1,b2)
    uu=zeros(m,n);
    uu(rd1,:)=b1;
    uu(:,rd2)=b2;
    rhs=M1*F*M2'-(K1*uu*M2'+M1*uu*K2');
    uu=V1(:,kd1)*((V1(kd1,kd1)'*rhs(kd1,kd2)*V2(kd2,kd2))./LL)*V2(:,kd2)';
end

y=y';
[xx,yy]=ndgrid(x,y);

F=ones(m,n);
b1=[0.2*sin(3*pi*y); 0*y];
b2=[(x<0).*sin(pi*x).^4, 0*x];
uu=greenF(F, 0*b1, 0*b2);

figure(1);
surf(xx,yy,uu);
colormap(jet(256));
shading interp;
camlight;
axis square manual; 
xlabel('x'); ylabel('y');
end


function [K, M, E, W, L, x]=GalerkinGLL(m,a,b)
% Stiffness and mass matrices for a Gauss-Legendre-Lobatto nodal basis over
% the interval [-1,1] with Robin BCs specified in a and b.
I=eye(m);                           % Identity
kd=2:m-1;                           % kept DOFs
rd=[1,m];                           % removed DOFs
[D,x,w]=legD(m);                    
D=D(end:-1:1,end:-1:1);             % Differentiation matrix
x=x(end:-1:1);                      % Collocation nodes
w=w(end:-1:1);                      % Quadrature weights
C=diag(a)*I(rd,:)+diag(b)*D(rd,:);	% Constraint operator

% Constrained basis
E=eye(m);
E(rd,kd)=-C(:,kd);
E(rd,:)=C(:,rd)\E(rd,:);

% Primed basis
DE=D*E;

% Mass matrix
V=VandermondeLeg(x);
Minv=(V*V');
M=E'*(Minv\E);

% Stiffness matrix
K=DE'*(diag(w)*DE)-E(rd,:)'*diag([1,-1])*DE(rd,:);

% Eigenfunctions
W=zeros(m);
[W(kd,kd),L]=eig(K(kd,kd), M(kd,kd), 'vector');
W(kd,kd)=normc(W(kd,kd), M(kd,kd));
W(rd,kd)=E(rd,kd)*W(kd,kd);
end
