function [] = HelmPolar(N,k)
[A1,A2,B1,r,phi]=chebLapPol(N,N);
A1=A1(2:end,2:end);
B1=B1(2:end,2:end);

xx=r*cos([0, phi]);
yy=r*sin([0, phi]);

figure(1);
[Vf,mu]=eigs(A2, k, 'sm');
Vf=bsxfun(@times, conj(Vf(1,:))./abs(Vf(1,:)), Vf);
mu=diag(mu);
I=eye(size(A1));
vv=zeros(N,N);

for j=1:k
    [Vr,lam]=eigs(A1+mu(j)*I, B1, k, 'sm');
    lam=sqrt(-diag(lam));
    Vr=bsxfun(@times, conj(Vr(end,:))./abs(Vr(end,:)), Vr);
    for i=1:k
        vv(2:end,:)=real(Vr(:,i)*Vf(:,j).');
        subplot(k,k,(i-1)*k+j);
        surf(xx,yy,[vv(:,end), vv]);
        axis square; axis off;
        shading interp;
        colormap(jet(256));
        text(-1,1.3,sprintf('%f', lam(i)));
        view(0,90);
    end
end
end