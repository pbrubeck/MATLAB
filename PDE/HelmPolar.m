function [] = HelmPolar(N,k)
[R2,Drr,Dff,r,phi]=chebLapPol(N,N);
Drr=Drr(2:end,2:end);
R2=R2(2:end,2:end);

xx=r*cos([0, phi]);
yy=r*sin([0, phi]);

figure(1);
[Vf,mu]=eigs(Dff, k, 'sm');
Vf=bsxfun(@times, conj(Vf(1,:))./abs(Vf(1,:)), Vf);
mu=diag(mu);
I=eye(size(Drr));
vv=zeros(N,N);

for j=1:k
    [Vr,lam]=eigs(Drr+mu(j)*I, R2, k, 'sm');
    lam=sqrt(-diag(lam));
    Vr=bsxfun(@times, conj(Vr(1,:))./abs(Vr(1,:)), Vr);
    for i=1:k
        vv(2:end,:)=real(Vr(:,i)*Vf(:,j).');
        subplot(k,k,(i-1)*k+j);
        surf(xx,yy,[vv(:,end), vv]);
        axis square; axis off;
        shading interp;
        colormap(gray(256));
        text(-1,1.3,sprintf('%f', lam(i)));
        view(0,90);
    end
end
end