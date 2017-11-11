function [] = wave1D(n)
% Solves the wave equation using Chebyshev collocation spectral method.

% Diff matrix, nodes and quadrature
[D,x]=chebD(n);
[w,~]=ClenshawCurtis(-1,1,n);
D2=D*D;

% Boundary conditions
a=[1,1]; % Dirichlet
b=[0,0]; % Neumann
kd=2:n-1;
rd=[1,n];
I=eye(n);
B=diag(a)*I(rd,:)+diag(b)*D(rd,:);
G=-B(:,rd)\B(:,kd);

% Schur complement
SA=D2(kd,kd)+D2(kd,rd)*G;

% Eigenmodes
S=zeros(n,n-2);
[S(kd,:), L]=eig(SA, 'vector');
S(rd,:)=G*S(kd,:);
omega=sqrt(-L);

% Initiall conditions
u=1-x.^2;
v=0*(1-x.^2);
bc=B*u;

% DC component
u0=zeros(n,1);
ub=G*u(kd)-u(rd);
u0(kd)=SA\(D2(kd,rd)*ub);
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
for i=1:nframes
    t=t+dt;
    u=u0+S*(cos(omega*t).*a0+(t*sinc(omega*t)).*b0);
    ut=  S*(-omega.*sin(omega*t).*a0+cos(omega*t).*b0);
    utt= S*(-omega.*(omega.*cos(omega*t).*a0+sin(omega*t).*b0));
    
    err=(utt-D*D*u);
    err(rd)=B*u-bc;
    set(h1, 'YData', u);
    set(h2, 'YData', err);
    drawnow;
    
    E(i)=(ut'*diag(w)*ut-(D*u)'*diag(w)*(D*u))/2;
end

figure(3);
plot(dt*(1:nframes), E-mean(E));
end

