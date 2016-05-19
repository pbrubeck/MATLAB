function [] = annulus(N)
% Dirichlet problem solved on the annulus
R1=2;
R2=4;
[D,x]=chebD(N);
r=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;

dt=2*pi/N;
phi=dt*(1:N);
xx=r*cos(phi);
yy=r*sin(phi);

% Laplacian
Drr=(diag(r)*D)^2;
d1=Drr(:,1); dn=Drr(:,end);
Dff=toeplitz([-pi^2/(3*dt^2)-1/6, 0.5*(-1).^(2:N)./sin(dt*(1:N-1)/2).^2]);
% Lap=kron(diag(r.^-2),eye(N))*(kron(Drr,eye(N))+kron(eye(N),Dff));

% Boundary conditions
g1=0*sin(10*phi);
g2=4*sin(5*phi);
F=zeros(N);
RHS=diag(r.^2)*F-dn*g1-d1*g2;

uu=zeros(N);
uu(N,:)=g1; uu(1,:)=g2;
uu(2:end-1,:)=sylvester(Drr(2:end-1,2:end-1), Dff, RHS(2:end-1,:));

err=diag(r.^-2)*(Drr*uu+uu*Dff)-F;
disp(mean(abs(err(:))));

figure(1); surfl(xx,yy,uu,'light'); 
shading interp; colormap(jet(256)); axis equal;
end