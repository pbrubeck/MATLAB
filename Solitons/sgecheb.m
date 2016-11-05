function [Q,x] = sgecheb(N, dt)
% Sine-Gordon Equation, Chebyshev linear propagator
% u_tt - u_xx = 0

% Domain substitution x=tan(th)
[D,th]=chebD(N); th=pi/2*th; D=diag(2/pi*cos(th).^2)*D;
x=tan(th);

% Linear propagator
[V,L]=eig(D(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N);
U(2:end-1,3:end)=V;
ND=null(D);
U(:,1:2)=ND(:,1:2);
Q1=real(U*diag(exp(dt*L))/U);
Q2=real(U*diag(dt*sinc(1i*dt*L))/U);
Q3=real(U*diag(exp(-dt*L))/U);

function [u,v]=sgechebprop(u,v)
    u=Q1*u+Q2*v;
    v=Q3*v;
end
Q=@sgechebprop;
end