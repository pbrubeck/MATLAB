function [] = domainDecomp( N, k )
% Schwarz Additive Method for the Helmholtz equation on the L-shaped
% membrane


[D,x]=chebD(N);
D2=D*D;
[x1,y1]=ndgrid((x-1)/2,x);
[x2,y2]=ndgrid(x,(x-1)/2);
g1=(x1==0) & (y1<0);
g2=(y2==0) & (x2<0);

function u=solvePoisson(F)
u1=zeros(N);
u2=zeros(N);
BC1=zeros(N);
BC2=zeros(N);
F1=reshape(F(1:N*N), [N,N]);
F2=reshape(F(1+N*N:2*N*N), [N,N]);


err=1; tol=1e-12;
while err>tol
    u2(g2)=interp2(x1.',y1.',u1.',x2(g2),y2(g2),'spline');
    BC2([1, N],:)=u2([1, N],:);
    BC2(:,[1, N])=u2(:,[1, N]);
    RHS2=F2-4*D2*BC2-BC2*D2';
    u2(2:N-1,2:N-1)=sylvester(4*D2(2:N-1,2:N-1), D2(2:N-1,2:N-1)', RHS2(2:N-1,2:N-1));
    
    u1(g1)=interp2(x2.',y2.',u2.',x1(g1),y1(g1),'spline');
    BC1([1, N],:)=u1([1, N],:);
    BC1(:,[1, N])=u1(:,[1, N]);
    RHS1=F1-D2*BC1-4*BC1*D2';
    u1(2:N-1,2:N-1)=sylvester(D2(2:N-1,2:N-1), 4*D2(2:N-1,2:N-1)', RHS1(2:N-1,2:N-1));
    
    err=norm(u2(g2)-interp2(x1.',y1.',u1.',x2(g2),y2(g2),'spline'));
end
u=[u1(:); u2(:)];
end

[V,D]=eigs(@solvePoisson,2*N*N,k,'sm');
v1=reshape(V(1:N*N,k), [N,N]);
v2=reshape(V(1+N*N:2*N*N,k), [N,N]);
lam=diag(D);

figure(2);
surf(x1,y1,v1); hold on;
surf(x2,y2,v2); hold off;
colormap(jet); shading interp; alpha(0.5);
axis square; view(2); 
title(sprintf('\\lambda_{%d}=%f',k,lam(k)));
end

