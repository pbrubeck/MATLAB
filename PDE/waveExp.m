function [] = waveExp(N)
% Solves the damped, forced wave equation
%
% u_tt + B u_t + A^2 u = F,
%
% where
% A^2 is the Laplace operator
% B = c.grad is the damping coefficient, advection operator
% F = F0*cos(q*t) is the inhomogenous forcing term
%
% given initial conditions
% u(x,y,0) = u0(x,y)
% u_t(x,y,0) = v0(x,y)
%
% subject to boundary conditions
% a*u + b*du/dn = constant


% Boundary conditions
a=[1,1;1,1]; 
b=[0,0;0,0]; 

% Advection
c=[0;0];

% Operators
rd=[1,N];
kd=2:N-1;
[D,x]=chebD(N); D2=D*D;
[xx,yy]=ndgrid(x);
[A1,G1]=setBC(D2, D, a(1,:), b(1,:));
[A2,G2]=setBC(D2, D, a(2,:), b(2,:));

% Initial conditions
% u0 = harmonic function + gaussian perturbation
% v0 = -c.grad(perturbation)
zz=xx+1i*yy;
hh=exp(-100*(xx.^2+yy.^2)/2);
U0=real(zz.^2)+hh;
V0=-(1*D*hh+0*hh*D');

% Eigenfunctions
W1=zeros(N); L1=zeros(N,1);
W2=zeros(N); L2=zeros(N,1);
W1(:,rd)=null(D2);
W2(:,rd)=null(D2);
[W1(kd,kd), L1(kd)] = eig(A1, 'vector');
[W2(kd,kd), L2(kd)] = eig(A2, 'vector');
W1(rd,kd)=G1*W1(kd,kd);
W2(rd,kd)=G2*W2(kd,kd);

% Eigenvalues, frequency and damping
[L1, L2] = ndgrid(L1, L2);
LL=L1+L2; 
BB=c(1)*sqrt(-L1)+c(2)*sqrt(-L2);
ww=sqrt(BB.^2-LL);
ww(rd,:)=0; ww(:,rd)=0;
BB(rd,:)=0; BB(:,rd)=0;

% Solution, metod of lines
UU=pinv(W1)*U0*pinv(W2)';
VV=(BB.*UU+pinv(W1)*V0*pinv(W2)')./ww;
VV(rd,:)=0; VV(:,rd)=0;
U=@(t) W1*(exp(-BB*t).*(UU.*cos(ww*t)+VV.*sin(ww*t)))*W2';

figure(1);
h=surf(xx,yy,U0);
colormap(jet(256));
zlim([-2,2]);
axis square manual;
camlight; shading interp; 

dt=2*pi/400;
nframes=1001;
for i=0:nframes-1
    uu=U(i*dt);
    set(h,'ZData',uu);
    drawnow;
end
end