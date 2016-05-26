function [] = annulus(N,M)
% Dirichlet problem solved on the annulus
R1=2;
R2=4;
[D,x]=chebD(N);
r=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;
Drr=(diag(r)*D)^2;
d1=Drr(:,1); dn=Drr(:,end);

[Dff, phi]=periodicD2(M);
xx=r*cos(phi);
yy=r*sin(phi);

% Laplacian
% Lap=kron(eye(M),diag(r.^-2))*(kron(eye(M),Drr)+kron(Dff,eye(N)));

% Boundary conditions
g1=0*cos(10*phi);
g2=4*sin(5*phi);
F=zeros(N,M);
RHS=diag(r.^2)*F-dn*g1-d1*g2;

uu=zeros(N,M);
uu(N,:)=g1; uu(1,:)=g2;
uu(2:end-1,:)=sylvester(Drr(2:end-1,2:end-1), Dff', RHS(2:end-1,:));

xx=[xx(:,end), xx];
yy=[yy(:,end), yy];
uu=[uu(:,end), uu];
figure(1); surfl(xx,yy,uu,'light'); 
shading interp; colormap(jet(256)); axis square;
end