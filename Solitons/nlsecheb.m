function [Q,x] = nlsecheb(N, dt)
% NonLinear Schrodinger Equation, Chebyshev linear propagator
% 1i u_t + 1/2 u_xx = 0

% Domain substitution x=tan(th), u=sec(th).*w
[D,th]=chebD(N); th=pi/2*th; D=2/pi*D;
x=tan(th);

T=-1/2*diag(cos(th).^4)*(eye(N)+D*D);
NT=null(T);
[V,L]=eig(T(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N);
U(:,1:2)=NT(:,1:2);
U(2:end-1,3:end)=V;
Q=U*(diag(exp(-1i*dt*L)))/U;
Q=Q.*(sec(th)*cos(th'));
end