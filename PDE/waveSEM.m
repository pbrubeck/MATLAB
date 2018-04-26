function [] = waveSEM(k,n)
% Solves the wave equation using Legedre collocation - weak Galerkin
% spectral element method and Pade approximation of the matrix exponential.

% Subscript 0 : element / subdomain
% Subscript 1 : whole grid
% Subscript 2 : phase space (scalar field + spacetime momentum) [u,h,p]

% Problem parameters
vel=0;
mass=0;        % mass=0 yields wave equation
L=4*pi;        % Spatial domain [-L,L]
T=4*L;         % Time domain [0,T]
dt=0.02;
nsteps=ceil(T/dt);
dt=T/nsteps;

N = (n-1)*k+1; % Total number of gridpoints
mask=(mod(0:(n-1)*(k-1),n-1)==0)';
% Element Diff matrix, nodes and quadrature
[D0,x0,w0]=legD(n);
V0=VandermondeLeg(x0);
M0=inv(V0*V0'); % Element mass
K0=diag(w0)*D0; % Element stiffness

% Element boundaries
xe=linspace(-L,L,k+1)';
J0=diff(xe)/2; % Jacobian
x1=kron(J0, x0)+kron((xe(1:end-1)+xe(2:end))/2, ones(n,1));
x1(n+1:n:end)=[];

% Spectral element assembly
J1=kron(J0,ones(n-1,1));
J1=J1(1:numel(mask));
M1=conv2(diag(mask.*J1), M0);
K1=conv2(diag(mask    ), K0);

% Degrees of freedom
rd=[N+1,3*N];         % Removed dofs
kd=setdiff(1:3*N,rd); % Kept dofs

% Constraint operator
B=zeros(2,3*N);
B(1,[N+1,2*N+1])=[1,1];
B(2,[2*N,3*N])=[1,-1];

% Schur complement basis, E1 satisfies homogenous BCs
E=eye(3*N,3*N);
E(rd,kd)=-B(:,kd);
E(rd,:)=B(:,rd)\E(rd,:);
E1=E(:,kd); %E2=E(:,rd);

% Time propagator
maskM1=[0,-1,0; mass,0,0; 0,0,0];
maskK1=[0,0,0; 0,0,-1; 0,-1,0];
Mfull=kron(eye(3),M1);
M2=E1'*Mfull*E1;
K2=E1'*(kron(maskM1,M1)+kron(maskK1,K1))*E1;
A=M2\K2;
U=expm(A*dt);

% [Lambda]=eig(U,'vector');
% figure(5);
% plot(Lambda,'.');

% Initial condition
u=-1/2*sin((-vel*x1-abs(x1))/(1-vel^2));
h= 1/2*cos((-vel*x1-abs(x1))/(1-vel^2)).*(1+vel*sign(x1))/(1-vel^2);
p= 1/2*cos((-vel*x1-abs(x1))/(1-vel^2)).*(vel+1*sign(x1))/(1-vel^2);
q=[u; h; p];

% Plot initial condition
figure(1); h1=plot(x1,[u,u]); xlim([-L,L]); ylim([-0.5,0.5]); title('Scalar field');
figure(2); h2=plot(x1,h); xlim([-L,L]); ylim([-0.5,0.5]); title('Enthalpy');
figure(3); h3=plot(x1,p); xlim([-L,L]); ylim([-0.5,0.5]); title('Momentum');

% Jump Force
xi=0;
Jumps=zeros(n+1,3);
Ju=zeros(n,1); Ju(2:4:end)=1;  Ju(4:4:end)=-1;
Jh=zeros(n,1); Jh(2:4:end)=1i; Jh(4:4:end)=-1i;
Jp=zeros(n,1); Jp(1:4:end)=1;  Jp(3:4:end)=-1;
Jumps(1:n,:)=[Ju,Jh,Jp];
diracM=[0,0,0];
diracK=[0,1,0];

JM=zeros(n,3);
JK=zeros(n,3);

JM(:,1)=-1i*Jumps(1:n,1)+Jumps(1:n,2);
JM(:,2)=-1i*Jumps(1:n,2);
JM(:,3)=-1i*Jumps(1:n,3);

JK(:,2)=Jumps(1:n,3);
JK(:,3)=Jumps(1:n,2);

f=jumpForce(xi,xe,x0,x1,M0,JM,diracM,-1)+jumpForce(xi,xe,x0,x1,K0,JK,diracK,0);
f=f(:)/2;

%Steady state
b=-(K2+1i*M2)\(E1'*f);

% Exact solution
uex=@(t) -1/2*sin((t-vel*x1-abs(x1-vel*t))/(1-vel^2));

% Timestepping
t=0;
for i=1:nsteps
    t=t+dt;
    q=E1*(U*real(q(kd)-b)+real(b*exp(-1i*dt)));
    b=b*exp(-1i*dt);
    
    set(h1,{'YData'},{q(1:N); uex(t)});
    set(h2, 'YData' , q(N+1:2*N));
    set(h3, 'YData' , q(2*N+1:3*N));
    drawnow;
end
end