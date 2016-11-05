function [] = nlsecheb(N, init)
% NLSE NonLinear Schrodinger Equation, Chebyshev method
% 1i u_t + 1/2 u_xx +|u|^2 u = 0
nframes=1024;
t0=-pi;
tf=pi;

dt=1/1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));

% Domain substitution x=tan(th), psi=sec(th).*u
[D,th]=chebD(N); th=pi/2*th; D=2/pi*D;
x=tan(th);

% Linear propagator
T=-1/2*diag(cos(th).^4)*(eye(N)+D*D);
[V,L]=eig(T(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N);
U(2:end-1,3:end)=V;
nH=null(T);
U(:,1:2)=nH(:,1:2);
Q=U*diag(exp(-1i*dt/2*L))/U;

nlse(init, Q, x, dt, t0, tf, sec(th));
end