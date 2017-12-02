function [] = waveDG(k,p)
% Solves the wave equation using Continuous Galerkin spectral method.
N=k*(p-1)+1;

xe=linspace(-1,1,k+1)'; % Element boundaries
[Dn,xn,wn]=legD(p); % Diff matrix, nodes and quadrature
J=diff(xe)/2; % Jacobian
x=kron(J, xn)+kron((xe(1:end-1)+xe(2:end))/2, ones(p,1));
x(p:p:end-1)=[];

% Mass matrix
V=VandermondeLeg(xn);
Mn=inv(V*V');
Mn=(Mn+Mn')/2;

% Stiffness matrix
Kn=Dn'*diag(wn)*Dn;
Kn=(Kn+Kn')/2;

% Galerkin Block-Matrices
K=zeros(N);
M=zeros(N);
D=zeros(N);
for i=1:k
    m=1+(i-1)*(p-1);
    K(m:m+p-1, m:m+p-1)=K(m:m+p-1, m:m+p-1)+Kn/J(i);
    M(m:m+p-1, m:m+p-1)=M(m:m+p-1, m:m+p-1)+Mn*J(i);
    D(m:m+p-1, m:m+p-1)=D(m:m+p-1, m:m+p-1)+Dn/J(i);
end

D(p:p-1:end-1,:)=D(p:p-1:end-1,:)/2;

% Boundary conditions
a=[1,1]; % Dirichlet
b=[0,0]; % Neumann
kd=2:N-1;
rd=[1,N];
I=eye(N);
B=diag(a)*I(rd,:)+diag(b)*D(rd,:);
G=-B(:,rd)\B(:,kd);

% Schur complement
SM=M(kd,kd)+M(kd,rd)*G;
SK=K(kd,kd)+K(kd,rd)*G;

% Eigenmodes
S=zeros(N,N-2);
[S(kd,:), L]=eig(SK, SM, 'vector');
S(rd,:)=G*S(kd,:);
omega=sqrt(L);

[~,id]=sort(omega);
figure(7);
plot(x,S(:,id(1:7)));

% Initiall conditions
u=1-x.^2;
v=0*(1-x.^2);
bc=B*u;

% DC component
u0=zeros(N,1);
ub=G*u(kd)-u(rd);
u0(kd)=SK\(K(kd,rd)*ub);
u0(rd)=G*u0(kd)-ub;

a0=S(kd,:)\(u(kd)-u0(kd));
b0=S(kd,:)\v(kd);

figure(1);
h1=plot(x,u);
axis manual;
ylim([-2,2]);

figure(2);
h2=plot(x,u);

t=0; tf=12;
nframes=1000;
dt=tf/nframes;
for i=1:nframes
    t=t+dt;
    u=u0+S*(cos(omega*t).*a0+(t*sinc(omega*t)).*b0);
    ut=  S*(-omega.*sin(omega*t).*a0+cos(omega*t).*b0);
    utt= S*(-omega.*(omega.*cos(omega*t).*a0+sin(omega*t).*b0));
    
    err=M*(utt-D*D*u);
    err(rd)=B*u-bc;
    set(h1, 'YData', D*u);
    set(h2, 'YData', err);
    drawnow;
end

end

