function [] = waveExp(N)
% Solves the damped, forced wave equation
%
% u_tt + B u_t + A^2 u = F,
%
% where
% A^2 is the Laplace operator
% B = damp+c.grad is the damping coefficient + advection operator
% F = F0*sin(q*t) is the inhomogenous forcing term
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
damp=1;

% Operators
rd=[1,N];
kd=2:N-1;
[D,x]=chebD(N); D2=D*D;
[xx,yy]=ndgrid(x);
[A1,G1]=setBC(D2, D, a(1,:), b(1,:));
[A2,G2]=setBC(D2, D, a(2,:), b(2,:));

% Initial conditions
% U0 = harmonic function
% V0 = p.grad(U0)
% F0 = forcing amplitude
x0=0.01;
w0=2*pi;
zz=xx+1i*yy;
hh=0*exp(-100*(xx.^2+yy.^2)/2);
U0=real(0*zz.^2)+hh;
V0=-(1*D*U0+0*U0*D');
F0=5000*(exp(-1000*((xx-x0).^2+yy.^2)/2)-exp(-1000*((xx+x0).^2+yy.^2)/2));

% Eigenfunctions
S1=zeros(N); L1=zeros(N,1);
S2=zeros(N); L2=zeros(N,1);
S1(:,rd)=null(D2);
S2(:,rd)=null(D2);
[S1(kd,kd), L1(kd)] = eig(A1, 'vector');
[S2(kd,kd), L2(kd)] = eig(A2, 'vector');
S1(rd,kd)=G1*S1(kd,kd);
S2(rd,kd)=G2*S2(kd,kd);

% Eigenvalues, frequency and damping
[L1, L2] = ndgrid(L1, L2);
AA=-(L1+L2);
BB=damp+c(1)*sqrt(-L1)+c(2)*sqrt(-L2);
WW=sqrt(AA+BB.^2);
WW(rd,:)=0; WW(:,rd)=0;
BB(rd,:)=0; BB(:,rd)=0;

% Solve force contribution F(:)=K\F0(:)
KK=(AA-w0^2).^2+4*w0^2*BB.^2;
FF=(pinv(S1)*F0*pinv(S2)')./KK;
FF(rd,:)=0; FF(:,rd)=0;

% Solution, method of lines
UU=pinv(S1)*U0*pinv(S2)';
VV=pinv(S1)*V0*pinv(S2)';

Y1=UU-2*w0*BB.*FF;
Y2=(VV-BB.*Y1-w0*(AA-w0^2).*FF)./WW;
Y2(rd,:)=0; Y2(:,rd)=0;
U=@(t) S1*(exp(-BB*t).*(Y1.*cos(WW*t)+Y2.*sin(WW*t))+...
          (sin(w0*t)*(AA-w0^2)+2*w0*cos(w0*t)*BB).*FF)*S2';

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