function [ ] = wavePade( n )
% Solves the wave equation using Legedre collocation - weak Galerkin
% spectral method and Pade approximation of the matrix exponential.

% Diff matrix, nodes and quadrature
[D,x,w]=legD(n);

% Boundary conditions
a=[1,1];  % Dirichlet
b=[1,-1]; % Neumann
kd=2:n-1; % kept DOFs
rd=[1,n]; % removed DOFs
I=eye(n);
B=diag(a)*I(rd,:)+diag(b)*D(rd,:); % Constraint operator
E=I;
E(rd,kd)=-B(:,kd);
E(rd,:)=B(:,rd)\E(rd,:);
E1=E(:,kd);  % Interior basis, Schur complement
E2=E(:,rd);  % Boundary basis
D1=D*E1;

% Mass matrix
V=VandermondeLeg(x);
Minv=(V*V');
M1=E1'*(Minv\E1);

% Stiffness matrix
K1=D1'*diag(w)*D1-E1(rd,:)'*diag([1,-1])*D1(rd,:);

% Force generator
F=(D1'*diag(w)*D(:,rd)-E1(rd,:)'*diag([1,-1])*D(rd,rd))/B(:,rd);

% Initial conditions
u=exp(-100*x.^2).*(1-100*x.^2);
v=0*(1-x.^2).^2;
bc=B*u;

% Force
f0=F*bc;
f1=[zeros(size(f0)); -M1\f0];

figure(1);
h1=plot(x,u); hold on;
h2=plot(x,u); hold off;
% xlim([-1,1]);
% ylim([-2,2]);

t=0; tf=2;
nsteps=1000; 
dt=tf/nsteps; dt^8

A = -M1\K1;
L=[zeros(size(A)), I(kd,kd); A, zeros(size(A))];
S=dt*L;

% Direct matrix exponentiation
U=expm(S);

% Implicit Pade for the first order reduction in time
II=eye(size(L));
HL=II+S*( 1/2*II+S*(3/28*II+S*( 1/84*II+S/1680)));
HR=II+S*(-1/2*II+S*(3/28*II+S*(-1/84*II+S/1680)));

% Schur complement of the block system
S=(dt^2)*A;
II=eye(size(A));
P=II+S*(3/28*II+S/1680);
R=dt*(1/2*II+S/84);
QL=P-R*(P\A*R);
QR=P+R*(P\A*R);
QM=2*A*R;

% Compute steady state
w=[u(kd); v(kd)];
w0=-L\f1;
w=w-w0;

v1 = zeros(size(v));
u0=E1*(-K1\f0)+E2*bc;
u=u-u0;

% % Eigenvalues of velocity update
% figure(2);
% [~,lambda]=eig(QR,QL,'vector');
% [~,id]=sort(lambda,'ascend');
% lambda=lambda(id);
% plot(lambda);
% 
% % Eigenvalues of full system
% figure(3);
% [~,z]=eig(HR,HL,'vector');
% [~,id]=sort(angle(z),'ascend');
% z=z(id);
% plot(z, '*');
% axis manual equal;
% xlim([-1,1]);
% ylim([-1,1]);


c=zeros(1,1);
% Timestepper
for i=1:nsteps
    t=t+dt;
    
    % Schur complement solver
    v1(kd) = QL\(QR*v(kd)+QM*u(kd));
    u(kd) = u(kd) + P\(R*(v1(kd)+v(kd)));
    v(kd) = v1(kd);
    u1 = u0 + E1*u(kd);
    c(1)=(v(kd)'*M1*v(kd)+u(kd)'*K1*u(kd))/2;
    set(h1, 'YData', 0*u1);
    
    % Direct solver
    w=HL\(HR*w);
    u1 = u1- (E1*(w(1:n-2)+w0(1:n-2))+E2*bc);
    set(h2, 'YData', u1);
    
    drawnow;
    %disp(c);
end

end