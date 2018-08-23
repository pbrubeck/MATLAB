function [x,w,D,Q] = nlsecheb(n, dt)
% NonLinear Schrodinger Equation, Chebyshev linear propagator
% 1i u_t + 1/2 u_xx = 0
rd=[1,n];
kd=setdiff(1:n,rd);
% Domain substitution x=tan(th), u=sec(th).*w
[Dt,th,w]=legD(n); th=pi/2*th; Dt=2/pi*Dt; w=pi/2*w;
x=tan(th);

K=-1/2*diag(cos(th).^4)*(eye(n)+Dt*Dt);
M=eye(n);
S=eye(n);
L=zeros(n,1);
[S(kd,kd),L(kd)]=eig(K(kd,kd),M(kd,kd),'vector');
Q=S*diag(exp(-1i*dt*L))/S;
Q=Q.*(sec(th)*cos(th'));
D=diag(cos(th).^2)*Dt;
w=(sec(th).^2).*w;
end