function [] = waveSEM(k,p)
% Solves the wave equation using Legedre collocation - weak Galerkin
% spectral element method and Pade approximation of the matrix exponential.

% Subscript 0 : element / subdomain
% Subscript 1 : whole grid

N = (p-1)*k+1; % Total number of gridpoints
mask=(mod(0:(p-1)*(k-1),p-1)==0)';

% Element Diff matrix, nodes and quadrature
[D0,x0,w0]=legD(p);
V0=VandermondeLeg(x0);
M0=inv(V0*V0'); % Element mass
K0=diag(w0)*D0; % Element stiffness

% Element boundaries
L=pi;
xe=linspace(-L,L,k+1)';
J0=diff(xe)/2; % Jacobian
x1=kron(J0, x0)+kron((xe(1:end-1)+xe(2:end))/2, ones(p,1));
x1(p+1:p:end)=[];

% Spectral element assembly
J1=kron(J0,ones(p-1,1));
J1=J1(1:numel(mask));
M1=conv2(diag(mask.*J1), M0);
K1=conv2(diag(mask    ), K0);

% Degrees of freedom classification
rd=[1,2*N];  % Removed
kd=2:2*N-1;  % Kept
kd1=1:N;     % Kept (Left)
kd2=N+(1:N); % Kept (Right)

% Mixed constraint operator
B=zeros(2,2*N);
B(1,[1,N+1])=[1,1];
B(2,[N,2*N])=[1,-1];

% Schur complement basis, E1 satisfies homogenous BCs
E=eye(2*N,2*N);
E(rd,kd)=-B(:,kd);
E(rd,:)=B(:,rd)\E(rd,:);
E1=E(:,kd); %E2=E(:,rd);

M2=E1'*kron([1,0; 0,1],M1)*E1;
K2=E1'*kron([0,1; 1,0],K1)*E1;

A=-M2\K2;

figure(5);
imagesc(log(abs(K2-K2')));
colormap(gray(256));
colorbar();


% Initial condition
u=1/2*cos(abs(x1));
v=1/2*cos(abs(x1)).*sign(x1);
w=[u;v];


% Force
xi=0;
jumps=zeros(p,1);
jumps(1:4:end)=-1;
jumps(2:4:end)=-1i;
jumps(3:4:end)= 1;
jumps(4:4:end)= 1i;
fu=jumpForce(xi,xe,x0,x1,K0,jumps);
fv=zeros(size(fu));

f=-N/(sqrt(L*3)*4)*[fu;fv];

% Transient
z=(A+1i*eye(size(A)))\(E1'*f);

figure(1);
h1=plot(x1,real(u)); hold on;
h2=plot(x1,real(v));
hold off;
xlim([-L,L]);
ylim([-0.5,0.5]);

figure(2);
h3=plot(x1,imag(u)); hold on;
h4=plot(x1,imag(v));
hold off;
xlim([-L,L]);
ylim([-0.5,0.5]);


figure(3);
h5=plot(x1,abs(real(u)./real(v)));
axis manual;


% Time propagator
dt=0.005;
U=expm(A*dt);
%[V,D]=eig(-K2, M2,'vector');


w0=E1'*w-z;
t=0;
T=6*L;
nsteps=ceil(T/dt);
dt=T/nsteps;
for i=1:nsteps
    t=t+dt;
    w0=U*w0;
    w(kd)=w0+exp(-1i*t)*z;    
    w1=E1*w(kd);
    
    set(h1,'YData',real(w1(kd1))); 
    set(h2,'YData',real(w1(kd2)));
    set(h3,'YData',imag(w1(kd1))); 
    set(h4,'YData',imag(w1(kd2)));
    set(h5,'YData',abs(real(w1(kd1))./real(w1(kd2))));
    drawnow;
end

end