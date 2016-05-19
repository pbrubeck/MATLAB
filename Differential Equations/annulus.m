function [] = annulus(N)
% Dirichlet problem solved on the annulus
R1=2;
R2=4;
[D,x]=chebD(N);
r=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;

phi=linspace(0,2*pi,N);
xx=r(:)*cos(phi);
yy=r(:)*sin(phi);

% Laplacian
Drr=(diag(r)*D)^2;
d1=Drr(:,1); dn=Drr(:,end);
Dff=(-2/N^2)*ifft(diag((1:N).^2))*dftmtx(N);
% Lap=kron(diag(r.^-2),eye(N))*(kron(Drr,eye(N))+kron(eye(N),Dff));

% Boundary conditions
g1=0*phi;
g2=4*sin(5*phi);
F=zeros(N);
RHS=diag(r.^-2)*F-dn*g1-d1*g2;

uu=zeros(N);
uu(N,:)=g1; uu(1,:)=g2;
uu(2:end-1,2:end-1)=sylvester(Drr(2:end-1,2:end-1),Dff(2:end-1,2:end-1).',RHS(2:end-1,2:end-1));

disp(norm(diag(r.^-2)*(Drr*uu+spPartialD(uu,2,2)),'fro'));

figure(1); surfl(xx,yy,real(uu),'light'); 
shading interp; colormap(jet(256))
end

