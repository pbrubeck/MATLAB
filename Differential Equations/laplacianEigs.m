function [] = laplacianEigs(N)
R1=1;
R2=2;
[D,x]=chebD(N);
r=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;

dt=2*pi/N;
phi=dt*(1:N);
xx=r*cos(phi);
yy=r*sin(phi);

% Laplacian
Drr=(diag(r)*D)^2;
Dff=toeplitz([-pi^2/(3*dt^2)-1/6, 0.5*(-1).^(2:N)./sin(dt*(1:N-1)/2).^2]);
Lap=kron(diag(r.^-2)*Drr,eye(N))+kron(diag(r.^-2),Dff);

[V,D]=eigs(Lap,2,'sm');
uu=reshape(V(:,end),[N,N]);

figure(1); surfl(xx,yy,uu,'light'); 
shading interp; colormap(jet(256));
end

