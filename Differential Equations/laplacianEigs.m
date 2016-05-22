function [] = laplacianEigs(N,k)
R1=0;
R2=1;
[D,x]=chebD(N);
r=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;
Drr=(diag(r)*D)^2;
[Dff,phi]=periodicD2(N);

xx=r*cos([0, phi]);
yy=r*sin([0, phi]);

Drr=Drr(2:end-1,2:end-1);
r=r(2:end-1);
% Laplacian (long method)
% Lap=kron(eye(N),diag(r.^-2)*Drr)+kron(Dff,diag(r.^-2));
% [V,d]=eigs(Lap,k,'sm');

figure(1);
[Vf,mu]=eigs(Dff, k, 'sm');
mu=diag(mu);
for j=1:k
    [Vr,lam]=eigs(Drr+mu(j)*eye(size(Drr)), diag(r.^2), k, 'sm');
    lam=diag(lam);
    for i=1:k
        vv=zeros(N);
        vv(2:end-1,:)=Vr(:,i)*Vf(:,j)';
        vv=[vv(:,end), vv];
        subplot(k,k,(i-1)*k+j);
        surfl(xx,yy,real(vv),'light'); 
        shading interp;
        colormap(jet(256));
        title(sprintf('%f', lam(i)));
        view(0,90);
    end
end
end