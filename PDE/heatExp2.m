function [] = heatExp2( N )
% Linear PDE with mixed boundary conditions
%
% u_t = u_xx + u_yy
% u + c*du/dn = constant
%
% Method of lines: discrete x, continous t
% u_t = D2*u -> u(t) = expm(t*D2)*u(0)
%
% D2 is singular, we discard first and last rows to introduce the boundary
% conditions and obtain the Schur complement of the new system. The
% exponential map of a SOLDO is computed through diagonalization: the
% eigenvectors are the union of those of the Schur complement and the
% nullspace of the SOLDO, and the eigenvalues are exp(dt*lambda).

c=0;
dt=0.0001;

[D,x]=chebD(N);
D2=D*D;
[A,G]=setBC(D2,D,[1,-1],c);

L=zeros(N,1);
V=zeros(N);
V(:,[1,N])=null(D2);
[V(2:N-1,2:N-1),L(2:N-1)]=eig(A,'vector');
V([1,N],2:N-1)=G*V(2:N-1,2:N-1);
Q=V*diag(exp(dt*L))*pinv(V);

[xx,yy]=ndgrid(x);
rr=hypot(xx,yy);
th=atan2(yy,xx);
uu=1*(rr<0.5 & abs(th)<0.9*pi);
uu([1,N],:)=G*uu(2:N-1,:);
uu(:,[1,N])=uu(:,2:N-1)*G';

figure(1);
h=surf(xx,yy,uu);
colormap(jet(256)); shading interp;
view(2); axis square manual;

nframes=1000;
for i=1:nframes
    uu=Q*uu*Q';
    set(h, 'ZData', uu);
    drawnow;
end
end