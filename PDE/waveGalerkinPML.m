function [ ] = waveGalerkinPML( n )
% Solves the wave equation using legedre collocation - weak Galerkin
% spectral method with a Perfectly Matched Layer.

% Diff matrix, nodes and quadrature
[D,x,w]=legD(n);

% Mass matrix
V=VandermondeLeg(x);
M=inv(V*V');
M=(M+M')/2;

% Stiffness matrix
K=D'*diag(w)*D;
K=(K+K')/2;

% Boundary conditions
a=[1,0]; % Dirichlet
b=[0,1]; % Neumann
kd=2:n-1;
rd=[1,n];
I=eye(n);
B=diag(a)*I(rd,:)+diag(b)*D(rd,:);
G=-B(:,rd)\B(:,kd);

% Schur complement
SM=M(kd,kd)+M(kd,rd)*G;
SK=K(kd,kd)+K(kd,rd)*G;

% Eigenmodes
S=zeros(n,n-2);
[S(kd,:), L]=eig(SK, SM, 'vector');
S(rd,:)=G*S(kd,:);
omega=sqrt(L);

% Initial conditions
u=3-x.^2+2*x;
v=0*(1-x.^2);
bc=B*u;

% DC component
u0=zeros(n,1);
ub=G*u(kd)-u(rd);
u0(kd)=SK\(K(kd,rd)*ub);
u0(rd)=G*u0(kd)-ub;

a0=S(kd,:)\(u(kd)-u0(kd));
b0=S(kd,:)\v(kd);

figure(1);
h1=plot(x,u);
axis manual;
ylim([-4,4]);

figure(2);
h2=plot(x,u);

t=0; tf=24;
nframes=200;
dt=tf/nframes;

E=zeros(nframes,1);
F=zeros(nframes,1);
for i=1:nframes
    t=t+dt;
    u=u0+S*(cos(omega*t).*a0+(t*sinc(omega*t)).*b0);
    ut=  S*(-omega.*sin(omega*t).*a0+cos(omega*t).*b0);
    utt= S*(-omega.*(omega.*cos(omega*t).*a0+sin(omega*t).*b0));
    
    err=M*(utt-D*D*u);
    err(rd)=B*u-bc;
    set(h1, 'YData', u);
    set(h2, 'YData', err);
    drawnow;
    
    E(i)=(ut'*M*ut+u'*K*u)/2;
end

figure(3);
plot(dt*(1:nframes), E-mean(E));
end