function [Q,x] = sgeherm(N, dt)
% Sine-Gordon Equation, Hermite linear propagator
% u_tt - u_xx = 0

[D,x]=hermD(N);
[V,L]=eig(D(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N);
U(2:end-1,3:end)=V;
U(:,1:2)=[0*x+1,x];
Q1=real(U*diag(exp(dt*L))/U);
Q2=real(U*diag(dt*sinc(1i*dt*L))/U);
Q3=real(U*diag(exp(-dt*L))/U);

function [u,v]=sgehermprop(u,v)
    u=Q1*u+Q2*v;
    v=Q3*v;
end
Q=@sgehermprop;
end