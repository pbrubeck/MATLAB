function []=discForce(n,xi)
% Degrees of freedom
kd=2:n-1;
rd=[1;n];

% Element Diff matrix, nodes and quadrature
[D0,x0,w0]=legD(n);
% Element Stiffness matrix
K0=D0'*diag(w0)*D0;

% Jacobians
h1=(1+xi)/2;
h2=(1-xi)/2;
% Centers
c1=(xi-1)/2;
c2=(xi+1)/2;
% Supergrid
xx=[c1+h1*x0; c2+h2*x0];
xx(n:n+1)=[xi-eps,xi+eps];
% Domain indices
d1=1:n;
d2=n+d1;

% Unitary jump in first derivative
nj=2;
jumps=zeros(nj,1);
jumps(2)=1;
[s1,s2]=piecewiseLagrange(x0,xi,jumps);
P=legC(x0,xx);

% Compute force
fdisc=-P'*[(1/h1)*K0*P(d1,:)*s1; (1/h2)*K0*P(d2,:)*s2];
f0=-P(n,:)';
f=f0+fdisc;

% Solve equation
u=zeros(n,1);
u(kd)=K0(kd,kd)\f(kd);

% Exact solution
a=(abs(1+xi)-abs(1-xi))/4;
b=1/2;
c=-0.25*((1-xi)*abs(1+xi)+(1+xi)*abs(1-xi));
uex=a*(x0-xi)+b*abs(x0-xi)+c;
err=abs(1-u./uex);
err(rd)=0;


% Interpolate jump function
jumps=zeros(nj,1);
jumps(4)=1;
Nq=2048;
xxx=linspace(-1,1,Nq)';
[s1,s2]=piecewiseLagrange(x0,xi,jumps);
P=legC(x0,xxx);
yy=[P(xxx<=xi,:)*s1; P(xxx>=xi,:)*s2];
xxx=[xxx(xxx<=xi); xxx(xxx>=xi)];

% Some plots
figure(1);
plot(xxx,yy,'b', x0,0*x0,'or');
title('Jump function');

figure(2);
plot(x0,fdisc);
title('Force');

figure(3);
plot(x0,u,'.', x0, uex);
title('Solution');

figure(4);
plot(x0, err);
title('Error');
end