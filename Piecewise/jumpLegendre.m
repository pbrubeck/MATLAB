function [] = jumpLegendre(n,m,xi)
% Legendre ODE with point source
[D0,x0]=legD(n);
VL=VandermondeGeg(1/2, x0);
VG=VandermondeGeg(3/2, x0);
M0=inv(VL*VL');
S0=VG\D0;
K0=S0'*S0;

[P,L]=eig(K0,M0,'vector');
A0=K0+L(m)*M0;


xe=[-1;1];
jumps=zeros(n,1);
jumps(1:2)=[0;1]/(1-xi^2);
for j=3:length(jumps)
    jumps(j)=((m+1)*m-(j-1)*(j-2))/(2*xi*(j-1))*jumps(j-1);
end
f=jumpForce(xi,xe,x0,x0,A0,jumps);


u0=A0\f;
u=u0+P(:,m);

figure(1);
plot(x0,u);
end

