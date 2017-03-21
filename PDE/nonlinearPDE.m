function [] = nonlinearPDE(N)
% Solves the nonlinear PDE using fixed point iteration
%
% Lu=f(u)
%
% where L is the Laplace operator and f(u) is a nonlinear function

% Sine-Gordon 2D steady state
% u_xx + u_yy = sin(u)

% Boundary conditions
a=[1,1;1,1]; 
b=[0,0;0,0]; 

% Operators
rd=[1,N];
kd=2:N-1;
[D,x]=chebD(N); D2=D*D;
[xx,yy]=ndgrid(x);
[A1,G1,H1,C1]=setBC(D2, D, a(1,:), b(1,:));
[A2,G2,H2,C2]=setBC(D2, D, a(2,:), b(2,:));

% Function
f=@(u) sin(u); 

% Initial guess
zz=xx+1i*yy;
uu=real(zz.^2);

% Newton-Raphson
i=0;
while(i<20)
    bc1=C1*uu;
    bc2=uu*C2';
    rhs=f(uu)-D2(:,rd)*uu(rd,:)-uu(:,rd)*D2(:,rd)';
    uu(kd,kd)=sylvester(A1, A2', rhs(kd,kd));
    uu(rd,:)=H1*bc1+G1*uu(kd,:);
    uu(:,rd)=bc2*H2'+uu(:,kd)*G2';
    i=i+1;
end

figure(1);
surf(xx,yy,uu);
colormap(jet(256));
shading interp;
camlight; axis square;

end