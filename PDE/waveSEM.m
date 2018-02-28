function [] = waveSEM(k,p)
% Solves the wave equation using Legedre collocation - weak Galerkin
% spectral element method and Pade approximation of the matrix exponential.

% Subscript 0 : element / subdomain
% Subscript 1 : whole grid

N = (p-1)*k+1; % Total number of gridpoints
mask=mod(0:(p-1)*(k-1),p-1)==0;

% Element Diff matrix, nodes and quadrature
[D0,x0,w0]=legD(p);
V0=VandermondeLeg(x0);
M0=inv(V0*V0'); % Element mass
K0=diag(w0)*D0; % Element stiffness
M0=(M0+M0')/2;
K0=(K0-K0')/2;

% Element boundaries
xb=linspace(-1,1,k+1)';
J=diff(xb)/2; % Jacobian
x1=kron(J, x0)+kron((xb(1:end-1)+xb(2:end))/2, ones(p,1));
x1(p+1:p:end)=[];

% Spectral element assembly
w1=conv(mask,w0);
M1=conv2(diag(mask), M0*J(1));
K1=conv2(diag(mask), K0);

% Boundary conditions
E1=zeros(2*N,2*N-2);
E1([1,N+1],1)=[1,-1*0];
E1(2:N-1,2:N-1)=eye(N-2);
E1(N+2:2*N-1,N:2*N-3)=eye(N-2);
E1([N,2*N],end)=[1*0,1];

Z1=zeros(size(M1));
M2=E1'*[M1, Z1; Z1, M1]*E1;
K2=E1'*[Z1, K1; K1, Z1]*E1;
A=-M2\K2;

% Initial condition
u=-100*x1.*exp(-100*(x1).^2/2);
v= 100*x1.*exp(-100*(x1).^2/2);
w=E1'*[u;v];

kd1=1:N;
kd2=N+(1:N);

figure(1);
h1=plot(x1,u); hold on;
h2=plot(x1,v);
hold off;
axis manual;

% Time propagation
dt=0.005;
U=expm(dt*A);


% % Eigenvalues of full system
% figure(3);
% [~,z]=eig(HR,HL,'vector');
% [~,id]=sort(angle(z),'ascend');
% z=z(id);
% plot(z, '*');
% axis manual equal;
% xlim([-1,1]);
% ylim([-1,1]);

T=20;
nsteps=ceil(T/dt);
for i=1:nsteps
    w=U*w;
        
    w1=E1*w;
    set(h1,'YData',w1(kd1)); 
    set(h2,'YData',w1(kd2));
    drawnow;
end

end

