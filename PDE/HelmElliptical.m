function [] = HelmElliptical(N,k)
f=1;
q=4; m=3; n=1024;
j=mod(m,2); r=(m-j+2)/2;
E(1:n-1)=q;
D=[j*(1+q), (j+2:2:2*n-2+j).^2];
if j==0
    E(1)=sqrt(2)*q;
end
[a,~]=trideigs(D,E);
a=a(r);

[Gu,Gv,Duu,Dvv,u,v]=chebLapEll(2*N+1,N);
Duu=Duu(2:end-1,2:end-1);
Gu=Gu(2:end-1,2:end-1);
u=u(1:N+1);

[Vu,lu]=eigs(Duu+a*eye(2*N-1),f^2*Gu,k,'sm');
[Vv,lv]=eigs(Dvv+a*eye(N),f^2*Gv,k,'sm');
lu=diag(lu);
lv=diag(lv);
idx=find(abs(lv+4*q)<1E-7);

Jem=zeros(N+1,k);
Jem(2:end,:)=Vu(1:N,:);
cem=bsxfun(@times, conj(Vv(1,:))./abs(Vv(1,:)), Vv);

figure(1); plot(u,(Jem));
figure(2); plot(v,cem(:,idx));

xx=f*cosh(u)*cos([0,v]);
yy=f*sinh(u)*sin([0,v]);

figure(3);
for i=1:k
    zz=(angle(1i*Jem(:,i)*cem(:,idx).'));
    subplot(4,k/4,i);
    surf(xx,yy,[zz(:,end), zz],'EdgeColor','none');
    colormap(gray(256));
    shading interp;
    axis square off;
    text(-f,1.5*f,sprintf('%f', lu(i)));
    view(0,90);
end

disp(lv(idx));
disp([lu,lv])

end