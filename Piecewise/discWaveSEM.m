function []=discWaveSEM(k,p,xi,omega0)
N = (p-1)*k+1; % Total number of gridpoints
mask=mod(0:(p-1)*(k-1),p-1)==0;

% Element Diff matrix, nodes and quadrature
[D0,x0,w0]=legD(p);
V0=VandermondeLeg(x0);
M0=inv(V0*V0'); % Element mass
K0=D0'*diag(w0)*D0; % Element stiffness
M0=(M0+M0')/2;
K0=(K0+K0')/2;

% Element boundaries
xe=linspace(-1,1,k+1)';
J=diff(xe)/2; % Jacobian
x1=kron(J, x0)+kron((xe(1:end-1)+xe(2:end))/2, ones(p,1));
x1(p+1:p:end)=[];

% Spectral element assembly
w1=conv(mask,w0);
M1=conv2(diag(mask), M0*J(1));
K1=conv2(diag(mask), K0);

% Degrees of freedom
kd=2:N-1;
rd=[1;N];

% Eigenmodes
V=zeros(N,N-2);
[V(kd,:),L]=eig(K1(kd,kd),M1(kd,kd),'vector');
freq=sqrt(L);

% Domain indices
d1=1:p;
d2=p+d1;

% Locate singularity
if any(xe==xi)
    [~,imax]=max(xe==xi);
    j=imax-1;
    jd=(1+(j-1)*(p-1)):(1+(j+1)*(p-1));
else
    [~,imax]=max(xe>xi);
    j=imax-1;
    jd=(1+(j-1)*(p-1)):(1+j*(p-1));
end

% Subdomain
xa=x1(jd(1));
xb=x1(jd(end));
% Jacobians
h0=(xb-xa)/2;
h1=(xi-xa)/2;
h2=(xb-xi)/2;
% Centers
c1=(xi+xa)/2;
c2=(xi+xb)/2;
% Supergrid
xs=[c1+h1*x0; c2+h2*x0];
% Canonical reference frame [-1,1]
xx=2/(xb-xa)*(xs-(xb+xa)/2);
xi0=2/(xb-xa)*(xi-(xb+xa)/2);

% Unitary jump in first derivative
nj=2;
jumps=zeros(nj,1);
jumps(1)=0;
jumps(2)=1;
[s1,s2]=piecewiseLagrange(x0,xi0,jumps);
P=legC(x0,xx);

% Compute force
fdisc=-P'*[(h0/h1)*K0*P(d1,:)*s1; (h0/h2)*K0*P(d2,:)*s2];
f0=fdisc-P(p+1,:)';
f=zeros(N,1);
f(jd)=h0*f0;

% Steady state
us=zeros(N,1);
us(kd)=(omega0^2*M1(kd,kd)+K1(kd,kd))\f(kd);

% Initial condition
ui=0*exp(-100*x1.^2/2);

ua=zeros(N-2,1);
ub=zeros(N-2,1);

ua=V(kd,:)\(ui(kd)-us(kd));

figure(1);
h=plot(x1, ui, 'k', 'LineWidth', 0.5);
ylim([-1,1]);
axis manual;


dt=0.01;
tf=100;
for t=0:dt:tf
    u=V*(cos(freq*t).*ua+sin(freq*t).*ub)+us*cos(omega0*t);
    set(h, 'YData', u);
    drawnow;
end


end