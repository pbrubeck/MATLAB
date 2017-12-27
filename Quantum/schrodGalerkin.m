function [] = schrodGalerkin(n)
% Solves the wave equation using Legedre collocation - weak Galerkin
% spectral method.

% Diff matrix, nodes and quadrature
[Dx,x,w]=legD(n);

% Boundary conditions
a=[1,1];  % Dirichlet
b=0*[-1,1]; % Neumann
kd=2:n-1;
rd=[1,n];
I=eye(n);
B=diag(a)*I(rd,:)+diag(b)*Dx(rd,:);
G=-B(:,rd)\B(:,kd);

% Schur complement
E=I(:,kd)+I(:,rd)*G;
SD=Dx(:,kd)+Dx(:,rd)*G;

% Mass matrix
V=VandermondeLeg(x);
Minv=(V*V');
SM=E'*(Minv\E);
SM=(SM+SM')/2;

% Stiffness matrix
SK=SD'*diag(w)*SD-G'*diag([1,-1])*SD(rd,:);
SK=(SK+SK')/2;

% Eigenmodes
S=zeros(n,n-2);
[S(kd,:), L]=eig(SK, SM, 'vector');
S(rd,:)=G*S(kd,:);
[L,id]=sort(L,'descend');
S=S(:,id);

% Force generator
F=-(G'*diag([1,-1])*Dx(rd,rd)+SD'*diag(w)*Dx(:,rd))/B(:,rd);

% Initiall conditions
sig=1/5;
k=3;
u=(pi*sig^2)^(-1/4)*exp(-(1/sig^2+1i*k^2)*(x.^2)/2);

% Normalization
u=u/sqrt(u'*(Minv\u));

% Force
bc=B*u;
f0=F*bc;

% Zero-energy component
u0=zeros(n,1);
u0(kd)=SK\f0;
u0(rd)=G*u0(kd)+B(:,rd)\bc;

% Eigenfunction expansion
a0=S(kd,:)\(u(kd)-u0(kd));

figure(1);
h1=plot(x,abs(u).^2);
ylim([0, 4]);drawnow; 

err=zeros(size(x));
figure(2);
h2=plot(x,err);

t=0; tf=10;
nframes=3000;
dt=tf/nframes;
for i=1:nframes
    t=t+dt;
    u=u0+S*(exp(-L/2i*t).*a0);
    ut=S*(-L/2i.*exp(-L/2i*t).*a0);
    
    err(kd)=2i*SM*ut(kd)+SK*(u(kd)-u0(kd));
    err(rd)=B*u-bc;
    
    set(h1, 'YData', abs(u).^2);
    set(h2, 'YData', abs(err));
    drawnow;
    H=real((u(kd)-u0(kd))'*SK*(u(kd)-u0(kd)))
    P=real((u(kd)-u0(kd))'*SM*(u(kd)-u0(kd)))
end
end

