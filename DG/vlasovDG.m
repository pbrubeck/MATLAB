function []=vlasovDG(k,p)

% One dimensional DG [0,L]
L=8;
xe=linspace(0,L,k+1)'; % Element boundaries
[Dx,xn,wn]=legD(p); % Diff matrix, nodes and quadrature
J=diff(xe)/2; % Jacobian
x=kron(J, xn)+kron((xe(1:end-1)+xe(2:end))/2,ones(p,1));

% Advection Numerical Flux
c=1; % velocity
s=1; % s=0 average, s=1 upwind
F=zeros(p,p+2);
F0=[1, sign(c)*s]*[0.5, 0.5; 0.5, -0.5];
F(1,1:2)=F0; F(end,end-1:end)=-F0;

% Mass matrix
V=VandermondeLeg(xn);
Minv=V*V';

% Stiffness matrix
K=zeros(p,p+2);
K(:,2:end-1)=Dx'*diag(wn);

% Stencil
S=Minv*(K+F);

% Galerkin Block-Matrix with periodic BCs
A=zeros(k*p);
for i=1:k
    m=1+(i-1)*p;
    A(m:m+p-1, 1+mod(m-2:m+p-1, k*p))=c*S/J(i);
end

% Force function
f=1./(x.^2+1e-1);
% Velocity
v=x;

% Operator splitting
[Z,d]=eig(A,'vector');
W=pinv(Z);
dt=0.005;
P1=exp(dt*d*v.');
P2=exp(dt*f*d.');
Q1=@(uu) real(Z*(P1.*(W*uu)));
Q2=@(uu) real((P2.*(uu*W.'))*Z.');

% Initial condition
[xx,vv]=ndgrid(x,v);
tau=1-eps;
bump=@(x) exp(1./(tau*x.^2-1));
uu=bump((2/L*xx-1)).*bump((2/L*vv-1));

% Plot
figure(1);
h0=surf(xx,vv,uu);
colormap(jet(256)); colorbar;
axis square; 
shading interp; 
view(2);%camlight;
drawnow;caxis manual;

T=5;
nframes=ceil(T/dt);
for i=1:nframes
    uu=Q2(Q1((uu)));
    set(h0,'ZData',uu);
    drawnow;
end
end

