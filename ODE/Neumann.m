function [] = Neumann(N)
% Solves u_xx=f(x) with free-fixed boundary conditions
[D,x]=chebD(N);
D2=D^2;
D2(end,:)=D(end,:);
d1=D2(:,1);
D2=D2(2:end,2:end);

% Boundary conditions
va=1;
ub=0;

f=exp(4*x);
f(end)=va;
rhs=f-d1*ub;
u=[ub; D2\rhs(2:end)];

figure(1);
plot(x,u);
end