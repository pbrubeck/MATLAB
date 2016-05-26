function [] = GutierrezVega(N,M)
% Solves Dirichlet BVP over an elliptical region
a=1;

R1=0.5;
R2=1;
[D,x]=chebD(N);
u=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;
Duu=D^2;
d1=Duu(:,1); dn=Duu(:,end);

[Dvv, v]=periodicD2(M);
xx=a*cosh(u)*cos(v);
yy=a*sinh(u)*sin(v);

% Boundary conditions
g1=0*sin(10*v);
g2=1*sin(5*v);
F=zeros(N,M);
RHS=a^2*(diag(sinh(u).^2)*F+F*diag(sin(v).^2))-dn*g1-d1*g2;

Psi=zeros(N,M);
Psi(N,:)=g1; Psi(1,:)=g2;
Psi(2:end-1,:)=sylvester(Duu(2:end-1,2:end-1), Dvv', RHS(2:end-1,:));

xx=[xx(:,end), xx];
yy=[yy(:,end), yy];
Psi=[Psi(:,end), Psi];
figure(1); surfl(xx,yy,Psi,'light'); 
shading interp; colormap(jet(256)); axis square;
end