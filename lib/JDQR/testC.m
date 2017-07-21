function []=testC(m,n,k)
rd1=[1,m];
rd2=[1,n];
kd1=2:m-1;
kd2=2:n-1;
kd=@(uu) reshape(uu(kd1,kd2),[],1);

E1=eye(m);
E2=eye(n);
[Dx,x]=chebD(m); Dxx=Dx*Dx;
[Dy,y]=chebD(n); Dyy=Dy*Dy; y=y';
[xx,yy]=ndgrid(x,y);

bc1=[1,0;1,0];
bc2=[1,0;1,0];
BC1=diag(bc1(:,1))*E1(rd1,:)+diag(bc1(:,2))*Dx(rd1,:);
BC2=diag(bc2(:,1))*E2(rd1,:)+diag(bc2(:,2))*Dy(rd2,:);
G1=-BC1(:,rd1)\BC1(:,kd1);
G2=-BC2(:,rd2)\BC2(:,kd2);

function uu=giveback(vv,kd1,kd2,rd1,rd2,G1,G2)
    uu=zeros(length(kd1)+length(rd1),length(kd2)+length(rd2));
    uu(kd1,kd2)=reshape(real(vv),length(kd1),length(kd2));
    uu(rd1,kd2)=G1*uu(kd1,kd2);
    uu(kd1,rd2)=uu(kd1,kd2)*G2';
    uu(rd1,rd2)=G1*uu(kd1,kd2)*G2';
end
gb=@(uu) giveback(uu,kd1,kd2,rd1,rd2,G1,G2);

A1=Dxx(kd1,kd1)+Dxx(kd1,rd1)*G1;
A2=Dyy(kd2,kd2)+Dyy(kd2,rd2)*G2;
[V1,L1]=eig(A1,'vector');
[V2,L2]=eig(A2,'vector');
[~,id]=sort(abs(L1),'ascend'); L1=L1(id); V1=V1(:,id);
[~,id]=sort(abs(L2),'ascend'); L2=L2(id); V2=V2(:,id);
LL=bsxfun(@plus, L1, L2.');

opA=@(uu) Dxx*uu+uu*Dyy';
afun=@(X) reshape(A1*reshape(X,size(V1,2),size(V2,2))+reshape(X,size(V1,2),size(V2,2))*A2.',[],1);
pfun=@(X) reshape(V1*((V1\reshape(X,size(V1,2),size(V2,2))/V2.')./LL)*V2.',[],1);
afunwrapper([],'set',afun,(m-2)*(n-2),pfun);

v0=kron(V1(:,1),V2(:,1));
figure(2);
[V,D]=jdqr('afunwrapper',numel(v0),k,0.0,struct('Disp',1,'LSolver','BiCGstab','Precond','afunwrapper','v0',v0));
lambda=real(diag(D));
[~,id]=sort(abs(lambda), 'ascend');
lambda=lambda(id);
V=V(:,id);

uu=gb(V(:,end));
uu=uu*(-1)^(max(uu(:))+min(uu(:))<0)/max(abs(uu(:)));

xq=linspace(-1,1,1024);
yq=linspace(-1,1,1024);
[xxx,yyy]=ndgrid(xq,yq);
uuu=interpcheb(interpcheb(uu,xq,1),yq,2);

figure(1);
surf(xxx,yyy,uuu);
colormap(jet(256)); colorbar;
shading interp; camlight;
axis square; view(2);
end