function [] = DirichletDisk(N,M)
% Dirichlet problem solved on the unit disk
[A1,A2,B1,r,phi]=chebLapPol(N,M);
xx=r*cos(phi);
yy=r*sin(phi);

% Boundary conditions
g=4*sin(5*phi);
F=zeros(N,M);
RHS=B1*F-A1(:,1)*g;

uu=zeros(N,M);
uu(1,:)=g;
uu(2:end,:)=sylvester(A1(2:end,2:end), A2', RHS(2:end,:));

xx=xx(:,[end 1:end]);
yy=yy(:,[end 1:end]);
uu=uu(:,[end 1:end]);
figure(1);
surf(xx,yy,uu); colormap(jet(256));
camlight; shading interp; axis square;
end

