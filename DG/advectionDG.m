function []=advectionDG(k,p)
xe=linspace(-1,1,k+1)'; % Element boundaries
[Dx,xn,wn]=legD(p); % Diff matrix, nodes and quadrature
J=diff(xe)/2; % Jacobian
x=kron(J, xn)+kron((xe(1:end-1)+xe(2:end))/2,ones(p,1));

% Advection Numerical Flux
c=1; % velocity
s=1; % s=0 average, s=1 upwind
F=zeros(p,p+2);
F0=[1, sign(c)*s]*[0.5, 0.5; 0.5, -0.5];
F(1,1:2)=F0; F(end,end-1:end)=-F0;

% Velocity function
c=reshape(1+4*(x+1)/2, p, k); c(:)=1;

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
    A(m:m+p-1, 1+mod(m-2:m+p-1, k*p))=diag(c(:,i))*S/J(i);
end

% Initial condition
u=exp(-100*(x+0.5).^2);

figure(1);
h1=plot(x,u);

% Time propagation
dt=0.005;
Q=expm(dt*A);

T=3;
nframes=ceil(T/dt);
for i=1:nframes
    u=Q*u;
    set(h1,'YData',u);
    drawnow;
end
end