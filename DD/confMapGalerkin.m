function [] = confMapGalerkin( m )
% Conformal mapping and Continuous Galerkin method
if(nargin==1)
    n=m;
end

kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];

% Vertices
v0=[1i;exp(1i*pi*7/6);exp(-1i*pi*1/6)];
%v0=2/(2-sqrt(2))*[1i;0;1];

% Sides
L=abs(v0([3,1,2])-v0([2,3,1]));
% Contact triangle
V=eye(3)+diag((sum(L)/2-L)./L([2,3,1]))*[-1,0,1; 1,-1,0; 0,1,-1];

z0=zeros(7,1);
z0([1,2,3])=v0;
z0([4,5,6])=V*v0;
z0(7)=(L'*v0)/sum(L); % Incenter

Z1=reshape(z0([7,4,5,1]),[2,2]);
Z2=reshape(z0([7,5,6,2]),[2,2]);
Z3=reshape(z0([7,6,4,3]),[2,2]);

% Grid
p=polygon(z0([1,5,7,4]));
f=rectmap(p, 1:4);
params=parameters(f);
xmin=min(real(params.prevertex));
xmax=max(real(params.prevertex));
ymin=min(imag(params.prevertex));
ymax=max(imag(params.prevertex));
dx=xmax-xmin;
dy=ymax-ymin;

% Boundary conditions
a=[1,1;1,1];
b=0*[0,1;0,1];

% Set differential operators
[K1,M1,E1,V1,L1,x]=GalerkinGLL(m,a(1,:),b(1,:));
[K2,M2,E2,V2,L2,y]=GalerkinGLL(n,a(2,:),b(2,:));
[Lx,Ly]=ndgrid(2/dx*L1,2/dy*L2); LL=Lx+Ly;
y=y';

% Jacobian determinant
[xx,yy]=ndgrid(xmin+dx*(x+1)/2, ymin+dy*(y+1)/2);
zz=xx+1i*yy;
J=dx/2*abs(evaldiff(f,zz)).^2;
J(rd1,rd2)=0;

function uu=greenF(F,b1,b2)
    uu=zeros(m,n);
    uu(rd1,:)=b1;
    uu(:,rd2)=b2;
    rhs=M1*(J.*F)*M2'-(2/dx*K1*uu*M2'+2/dy*M1*uu*K2');
    uu(kd1,kd2)=V1(kd1,kd1)*((V1(kd1,kd1)'*rhs(kd1,kd2)*V2(kd2,kd2))./LL)*V2(kd2,kd2)';
    uu=E1*uu*E2';
end

F=ones(m,n);
b1=[0.2*sin(3*pi*y); 0*y];
b2=[(x<0).*sin(pi*x).^4, 0*x];
uu=greenF(F, 0*b1, 0*b2);

% Conformal mapping ww=f(zz)
b1=f(zz([1 end],:)); b2=f(zz(:,[1 end]));
ww=greenF(zeros(m,n), b1, b2);

% Plot solution
figure(1);
surf(real(ww),imag(ww),uu); colormap(jet(256));
shading interp; camlight; view(2); 
xlabel('x'); ylabel('y'); axis square manual;
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
