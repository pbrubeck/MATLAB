function [] = DirichletRing(N,M)
% Dirichlet problem solved on the annulus
R1=2;
R2=4;
[D,x]=chebD(N);
r=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;
Drr=(diag(r)*D)^2;

[Dff, phi]=periodicD2(M);
xx=r*cos(phi);
yy=r*sin(phi);

% Boundary conditions
b=[4*sin(5*phi); 0*cos(10*phi)];
F=zeros(N,M);
RHS=diag(r.^2)*F-Drr(:,[1 end])*b;

uu=zeros(N,M);
uu([1 end],:)=b;
uu(2:end-1,:)=sylvester(Drr(2:end-1,2:end-1), Dff', RHS(2:end-1,:));

xx=xx(:,[end 1:end]);
yy=yy(:,[end 1:end]);
uu=uu(:,[end 1:end]);
figure(1); surfl(xx,yy,uu,'light'); 
shading interp; colormap(jet(256)); axis square;
end