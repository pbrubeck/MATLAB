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
dt=0.001;

rd=[1,N];
kd=2:N-1;
[D,x]=chebD(N);
D2=D*D;

I=eye(N);
B=diag([1,-1])*I([1,end],:)+diag([c,c])*D([1,end],:); 
[A,G]=setBC(D2,D,[1,-1],c);

L=zeros(N,1);
V=zeros(N);
V(:,rd)=null(D2);
[V(kd,kd),L(kd)]=eig(A,'vector');
V(rd,kd)=G*V(kd,kd);

[L1,L2]=ndgrid(L);
LL=L1+L2;
LL(L1.*L2==0)=0; % what a move !!
Q=exp(dt*LL);

[xx,yy]=ndgrid(x);
rr=hypot(xx,yy);
th=atan2(yy,xx);
uu=1*(rr<0.5 & abs(th)<0.9*pi);

sig=0.5;
uu=uu+exp(-((xx-0.4).^2+yy.^2)/(2*sig^2)+1i*xx);

figure(1);
h=surf(xx,yy,abs(uu).^2);
colormap(jet(256)); camlight; shading interp;
axis square manual;

nframes=1000;
for i=1:nframes
    uu=V*(Q.*(V\uu/V'))*V';
    set(h, 'ZData', abs(uu).^2);
    drawnow;
end
end