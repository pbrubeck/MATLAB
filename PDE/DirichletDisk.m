function [] = DirichletDisk(N,M)
% Dirichlet problem solved on the unit disk
[R2,Drr,Dff,r,phi]=chebLapPol(N,M);
d1=Drr(:,1);
xx=r*cos(phi);
yy=r*sin(phi);

% Boundary conditions
g=4*sin(5*phi);
F=zeros(N,M);
RHS=R2*F-d1*g;

uu=zeros(N,M);
uu(1,:)=g;
uu(2:end,:)=sylvester(Drr(2:end,2:end), Dff', RHS(2:end,:));

xx=[xx(:,end), xx];
yy=[yy(:,end), yy];
uu=[uu(:,end), uu];
figure(1); surfl(xx,yy,uu,'light'); 
shading interp; colormap(jet(256)); axis square;
end

