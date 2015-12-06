function [ uu ] = LaplaceBVP(N)
[D,x]=chebD(N);
[xx, yy]=meshgrid(x); x=xx(:); y=yy(:);
D2=D^2; I=eye(N);
L=kron(D2, I)+kron(I, D2);

b=find(abs(x)==1 | abs(y)==1);
L(b,:)=zeros(4*(N-1), N^2);
L(b,b)=eye(4*(N-1));
rhs=zeros(N^2, 1);
rhs(b)=(y(b)==1 & x(b)<0).*sin(pi*x(b)).^4+0.2*(xx(b)==1).*sin(3*pi*yy(b));

tic; u=L\rhs; uu=reshape(u, N, N); toc
figure(1); surf(xx,yy,uu); shading interp;
colormap(jet(256)); whitebg('k');
end

