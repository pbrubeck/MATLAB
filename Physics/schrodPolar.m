function [] = schrodPolar(N, k)
R1=1;
R2=3;
[D,x]=chebD(N);
r=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;
Drr=(diag(r)*D)^2;
[Dtt,th]=fourD2(N);
th=[0,th];
xx=r*cos(th);
yy=r*sin(th);

Drr=Drr(2:end-1,2:end-1);
r=r(2:end-1);

V=200*((r-R1)/(R2-R1)).^2;
figure(2); plot(r,V);

[Phi,mu]=eigs(Dtt, k, 'sm');
mu=diag(mu);
Phi=Phi/sqrt(2*pi/N);

[x,w]=ClenshawCurtis(R1,R2,N);
W=diag(w(2:end-1));

E=zeros(k,k);
RR=zeros(N,k,k);
for j=1:k
    [R,lam]=eigs(Drr+mu(j)*eye(size(Drr))-diag(V.*r.^2), diag(r.^2), k, 'sm');
    E(:,j)=-diag(lam); 
    RR(2:end-1,:,j)=R*diag(1./sqrt(r'*W*(R.^2)));
end
c=(rand(size(E))+1i*rand(size(E)))./E;
c=c/norm(c,'fro');

frames=5000;
tmax=0.2;
dt=1/frames;
Psi=zeros(N,N);
zz=zeros(N,N+1);

h=surf(xx,yy,zz); view(0,90);
shading interp; colormap(jet(256));
zlim([0, 1]);
for t=0:dt:tmax
    Psi(:,:)=0;
    a=c.*exp(-1i*E*t);
    for j=1:k
        Psi=Psi+(RR(:,:,j)*a(:,j))*Phi(:,j)';
    end
    zz=abs(Psi).^2;
    set(h,'ZData',[zz(:,end),zz]);
    drawnow;
end
end