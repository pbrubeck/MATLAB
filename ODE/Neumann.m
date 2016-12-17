function [] = Neumann(N,k)
% Example on Neumann boundary conditions
[~,w]=ClenshawCurtis(-1,1,N);
[D,x]=chebD(N);
D2=D*D;
[A,G,H]=setBC(D2,D,[0,1],[1,0]);

% Helmholtz equation
U=zeros(N,N-2);
[U(2:N-1,:),L]=eig(A,'vector');
[L,id]=sort(L,'descend');
U([1,N],:)=G*U(2:N-1,:);
U=normc(U, w);

figure(1);
plot(x,U(:,id(k)));
disp(L(k)*4/(pi^2));

% Poisson Equation
f=exp(-4*x(2:N-1));
bc=[0;1];

u=zeros(N,1);
u([1,N])=H*bc;
u(2:N-1)=A\(f-D2(2:N-1,[1,N])*u([1,N]));
u([1,N])=u([1,N])+G*u(2:N-1);

figure(2);
plot(x, u);
end