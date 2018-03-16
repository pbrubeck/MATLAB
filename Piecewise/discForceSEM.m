function []=discForceSEM(k,p,xi)
N = (p-1)*k+1; % Total number of gridpoints
mask=(mod(0:(p-1)*(k-1),p-1)==0)';

% Element Diff matrix, nodes and quadrature
[D0,x0,w0]=legD(p);
V0=VandermondeLeg(x0);
M0=inv(V0*V0');
K0=D0'*diag(w0)*D0;

% Element boundaries
xe=linspace(-1,1,k+1)';
J0=diff(xe)/2; % Jacobian
x1=kron(J0, x0)+kron((xe(1:end-1)+xe(2:end))/2, ones(p,1));
x1(p+1:p:end)=[];

% Spectral element assembly
J1=kron(J0,ones(p-1,1));
J1=J1(1:numel(mask));
M1=conv2(diag(mask.*J1), M0);
K1=conv2(diag(mask./J1), K0);

% Right hand side
f=jumpForce(xi,xe,x0,x1,K0,[1;-1]);

% Degrees of freedom
kd=2:N-1;
rd=[1;N];

% Solve equation
u=zeros(N,1);
u(kd)=K1(kd,kd)\f(kd);

% Exact solution
a=(abs(1-xi)-abs(1+xi))/4;
b=-1/2;
c=0.25*((1-xi)*abs(1+xi)+(1+xi)*abs(1-xi));
uex=a*(x1-xi)+b*abs(x1-xi)+c;
err=abs(1-u./uex);
err(rd)=nan;

% Interpolate jump function
jumps=[1;1];
Nq=2049;
xxx=linspace(-1,1,Nq)';
[s1,s2]=piecewiseLagrange(x0,xi,jumps);
C=legC(x0,xxx);
yy=[C(xxx<=xi,:)*s1; C(xxx>=xi,:)*s2];
xxx=[xxx(xxx<=xi); xxx(xxx>=xi)];

% Some plots
figure(1);
plot(xxx,yy,'b', x0,0*x0,'or');
title('Jump function');

figure(2);
plot(x1,f);
title('Force');

figure(3);
plot(x1,u,'.r', x1, uex,'k');
title('Solution');

figure(4);
plot(x1, err);
title('Error');
end