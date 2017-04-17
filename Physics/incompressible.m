function [] = incompressible(n)
% Navier Stokes equataion for incompresible flow, Chebyshev spectral method

% Simulation parameters
nu=0.01;   % Viscosity
dt=0.01/n^2; % Timestep
nframes=100000;

% Differential operators
rd=[1,n];
kd=2:n-1;
[D,x]=chebD(n); D2=D*D;
[xx,yy]=ndgrid(x); y=x';

% Poisson solver, imposition of boundary conditions
B=D(rd,:);
G=-B(:,rd)\B(:,kd);
A=D2(kd,kd)+D2(kd,rd)*G;
S=zeros(n,n-2);
[S(kd,:),L]=eig(A,'vector');
S(rd,:)=G*S(kd,:);
W=inv(S(kd,:));
[L1,L2]=ndgrid(L); LL=L1+L2;
LL(abs(LL)<1e-5)=inf;

V=zeros(n,n-2); 
V(kd,:)=eye(n-2);
V(rd,:)=G;
N=zeros(n,2); 
N(rd,:)=inv(B(:,rd));
P=V/(V'*V)*V';
Q=eye(n)-P;

% Solenoidal projection operator
function [u]=solenoidal(u)
    b1=real(u(rd,:));
    b2=imag(u(:,rd));
    F=b1*Q'*N+N'*Q*b2;
    phi0=N*sylvester(N'*Q*N, N'*Q*N, F)*N';
    phib=phi0+(N*b1-phi0)*P+P*(b2*N'-phi0);
    div=D*real(u)+imag(u)*D';
    rhs=div-D2*phib-phib*D2';
    phi=phib+S*((W*rhs(kd,kd)*W')./LL)*S';
    u=u-D*phi-1i*(phi*D');
end

% Partial time derivative
function [v]=partialTime(u)
    adv=real(u).*(D*u)+imag(u).*(u*D');
    dif=D2*u+u*D2';
    v=solenoidal(-adv+nu*dif);
end

% Runge-Kutta
function [u]=solveRK4(u)
    k1=dt*partialTime(u     );
    k2=dt*partialTime(u+k1/2);
    k3=dt*partialTime(u+k2/2);
    k4=dt*partialTime(u+k3  );
    u=u+(k1+2*k2+2*k3+k4)/6;
end

% Initial condition
zz=xx+1i*yy;
u=-2000i*(xx+2i*(yy-1/2)).*exp(-3*abs(zz+1i/2).^2);
u=solenoidal(u);

[xq,yq]=ndgrid(linspace(-1,1,50));
uq=interp2(yy,xx,u,yq,xq,'spline');
uq=uq./abs(uq);

figure(1); clf; hold on;
h1=imagesc(x,y,abs(u));
colormap(jet(256)); colorbar;
axis square manual; 
xlabel('x'); ylabel('y');
xlim([-1,1]); ylim([-1,1]);
h2=quiver(xq,yq,real(uq),imag(uq),'w');
hold off;

% Euler method
for i=1:nframes
    u=solveRK4(u);
    u(rd,:)=1i*imag(u(rd,:));
    u(:,rd)=real(u(:,rd));
    u(rd,rd)=0;
    u=solenoidal(u);
    
    uq=interp2(yy,xx,u,yq,xq,'spline');
    uq=uq./abs(uq);
    set(h1,'CData',abs(u));
    set(h2,'UData',real(uq));
    set(h2,'VData',imag(uq)); 
    drawnow;
end
end