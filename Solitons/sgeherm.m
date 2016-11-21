function [Q,x] = sgeherm(N, dt)
% Sine-Gordon Equation, Hermite linear propagator
% u_tt - u_xx = 0

[D,x]=hermD(N);
[V,L]=eig(D(2:end-1,2:end-1),'vector');
L=[0;0;L];
S=zeros(N);
S(2:end-1,3:end)=V;
S(:,1:2)=[0*x+1,x];
Q1=real(S*diag(exp(dt*L))/S);
Q2=real(S*diag(dt*sinc(1i*dt*L))/S);
Q3=real(S*diag(exp(-dt*L))/S);

function [u,v]=sgehermprop(u,v)
    u=Q1*u+Q2*v;
    v=Q3*v;
end
Q=@sgehermprop;
end