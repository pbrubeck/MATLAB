function [x,w,D,Q] = nlsecheb(N, dt)
% NonLinear Schrodinger Equation, Chebyshev linear propagator
% 1i u_t + 1/2 u_xx = 0

% Domain substitution x=tan(th), u=sec(th).*w
[Dt,th]=chebD(N); th=pi/2*th; Dt=2/pi*Dt;
x=tan(th);
[~,w]=ClenshawCurtis(-pi/2,pi/2,N);
w=(sec(th).^2).*w(:);

T=-1/2*diag(cos(th).^4)*(eye(N)+Dt*Dt);
NT=null(T);
[V,L]=eig(T(2:end-1,2:end-1),'vector');
L=[0;0;L];
S=zeros(N);
S(:,1:2)=NT(:,1:2);
S(2:end-1,3:end)=V;
Q=S*(diag(exp(-1i*dt*L)))/S;
Q=Q.*(sec(th)*cos(th'));

D=diag(cos(th).^2)*Dt;
end