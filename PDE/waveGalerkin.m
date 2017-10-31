function [ ] = waveGalerkin( n )
% Solves the wave equation using legedre collocation - weak Galerkin
% spectral method.

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
a=[1,1]; % Dirichlet
b=[0,0]; % Neumann
kd=2:n-1;
rd=[1,n];
I=eye(n);
B=diag(a)*I(rd,:)+diag(b)*D(rd,:);
H=inv(B(:,rd));
N=zeros(n,2); N(rd,:)=H;
G=-B(:,rd)\B(:,kd);

% Schur complement
SM=M(kd,kd)+M(kd,rd)*G;
SK=K(kd,kd)+K(kd,rd)*G;

% Eigenmodes
[S, L]=eig(SK, SM, 'vector');
omega=sqrt(L);

% Initiall conditions
u=(exp(-100*x.^2/2));
v=0*u;

% DC component
u0=zeros(n,1);
u0(rd)=u(rd)-G*u(kd);
u0(kd)=-SK\(K(kd,rd)*u0(rd));

v=v-N*B*v;
a0=S\(u(kd)-u0(kd));
b0=S\v(kd);

ut=zeros(n,1);
utt=zeros(n,1);

figure(1);
h1=plot(x,u);
axis manual;
ylim([-4,4]);

figure(2);
h2=plot(x(kd),u(kd));

t=0; tf=12;
nframes=128;
dt=tf/nframes;

P=zeros(n-2,n);
P(:,kd)=eye(n-2);
E=zeros(nframes,1);
c=zeros(nframes,1);
for i=1:nframes
    t=t+dt;
    u(kd)=u0(kd)+S*(cos(omega*t).*a0+(t*sinc(omega*t)).*b0);
    u(rd)=u0(rd)+G*u(kd);
    ut(kd)=S*(-omega.*sin(omega*t).*a0+cos(omega*t).*b0);
    ut(rd)=G*ut(kd);
    utt(kd)=S*(omega.*(-omega.*cos(omega*t).*a0+sin(omega*t).*b0));
    utt(rd)=G*utt(kd);
    
    set(h1, 'YData', u);
    set(h2, 'YData', P*(utt-D*D*u));
    drawnow;
    
    Eb=diag(abs(a./b))*(u(rd)); %sum(Eb(~isnan(Eb)))
    E(i)=(ut'*M*ut+u'*K*u);
end

figure(3);
plot(1:nframes, E);
end

