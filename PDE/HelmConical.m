function [] = HelmConical(N, k)
a=2;
b=sqrt(2);

[D,x]=chebD(N);

[Dr, r]=chebD(2*N);
Lr=Dr*diag(r.^2)*Dr;
Lr=Lr(1:N,1:N)+Lr(1:N,end:-1:N+1);
r=r(1:N);

u=(a-b)/2*(x+1);
Du=2/(a-b)*D;
Lu=diag(-2*u.^3+(a^2+b^2)*u)*Du+diag((u.^2-a^2).*(u.^2-b^2))*Du^2;

v=b*x;
Dv=1/b*D;
Lv=diag(2*v.^3-(a^2+b^2)*v)*Dv+diag((v.^2-a^2).*(v.^2-b^2))*Dv^2;

[Vr, m]=eigs(Lr(2:end,2:end), k, 'sm');
R=zeros(N,k);
R(2:end,:)=Vr;
for j=1:k
    [Vu, lu]=eigs(Lu(2:end-1,2:end-1)+m(j)*diag(u(2:end-1).^2), k, 'sm');
    U=zeros(N,k);
    U(2:end-1,:)=Vu;
    [Vv, lv]=eigs(Lv(2:end-1,2:end-1)-m(j)*diag(v(2:end-1).^2), k, 'sm');
    V=zeros(N,k);
    V(2:end-1,:)=Vv;
end

figure(1); plot(r, R);
figure(2); plot(u, U);
figure(3); plot(v, V);
lu=diag(lu);
lv=diag(lv);
disp([lu, lv]);
end

