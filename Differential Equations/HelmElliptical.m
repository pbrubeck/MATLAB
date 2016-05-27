function [] = HelmElliptical(N,k)
[F,G,Duu,Dvv,u,v]=chebLapEll(N,N);
Duu=Duu(2:end,2:end);
F=F(2:end,2:end);
Z=zeros(N-1,N);

a=1.001;
A=[Duu-a*eye(N-1), Z; Z', Dvv+a*eye(N)];
B=[F, Z; Z', G];

[V,lam]=eigs(A,B,k,'sm');
lam=diag(lam);
Jem=zeros(N,k);
Jem(2:end,:)=V(1:N-1,:);
cem=V(N:end,:);
figure(1); plot(u,real(Jem));
figure(2); plot(v,real(cem));

xx=cosh(u)*cos([0,v]);
yy=sinh(u)*sin([0,v]);
%surf(xx,yy,zz);

figure(3);
for i=1:k
    if(norm(Jem(:,i))<1E-10)
        Jem(:,i)=0;
    end
    zz=real(Jem(:,i)*cem(:,i).');
    subplot(4,k/4,i);
    surf(xx,yy,[zz(:,end), zz],'EdgeColor','none');
    axis square; axis off;
    colormap(jet(256));
    text(-1,1.3,sprintf('%f', lam(i)));
    view(0,90);
end
end