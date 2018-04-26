function [erru,errv,errp,h,pcalls] = femstokes(gmshfile,n)
restart=20;
tol=1E-6;
maxit=10;

% Canonical domain
z0=[0; 1; 1i; 0.5; 0.5+0.5i; 0.5i];
x0=real(z0);
y0=imag(z0);
P0=[ones(size(x0)), 2*x0, 2*y0, x0.^2, 2*x0.*y0, y0.^2]; 
R=full(sparse(1:9,[1,2,3,2,4,5,3,5,6],1,9,6));
% Quadrature in canonical domain
[zq,wq]=trigauss(n); 
% Vandermonde matrix on quadrature nodes
Vq=[ones(size(zq)), real(zq), imag(zq)];
% Interpolation to quadrature points
Phi1=Vq*[1,0,0;-1,1,0;-1,0,1];
Phi2=(kron(Vq,ones(1,3)).*kron(ones(1,3),Vq))*(R/P0);
% Derivative of cardinal basis functions, evaluate linear
Dx_Phi2=Vq*(full(sparse(1:3,[2,4,5],2,3,6))/P0);
Dy_Phi2=Vq*(full(sparse(1:3,[3,5,6],2,3,6))/P0);

[tri, vert, bnd] = loadgmsh(gmshfile);
T=size(tri,1);
N2=size(vert,1);

id1=unique(reshape(tri(:,1:3),[],1));
big2small=zeros(N2,1);
big2small(id1)=1:numel(id1);
N1=size(id1,1);

m1=3;
ei1=zeros(m1*m1,T);
ej1=zeros(m1*m1,T);
M1=zeros(m1*m1,T);
K1=zeros(m1*m1,T);

m2=6;
ei2=zeros(m2*m2,T);
ej2=zeros(m2*m2,T);
M2=zeros(m2*m2,T);
K2=zeros(m2*m2,T);

ei3=zeros(m2*m1,T);
ej3=zeros(m2*m1,T);
DX=zeros(m2*m1,T);
DY=zeros(m2*m1,T);

h=0;
for e=1:T
    nodes=tri(e,:);
    Pe=[ones(3,1), real(vert(nodes(1:3))), imag(vert(nodes(1:3)))];
    ze=vert(nodes);
    Area=abs(det(Pe))/2;
    h=max(h, sqrt(2*Area));
    
    % Linear elements (pressure)
    C=inv(Pe);
    grad=C(2:3,:);
    Ke=Area*(grad'*grad);
    Me=Area/12*(ones(3)+eye(3));
    [ni,nj]=ndgrid(big2small(nodes(1:m1)));
    ei1(:,e)=ni(:);
    ej1(:,e)=nj(:);
    M1(:,e)=Me(:);
    K1(:,e)=Ke(:);    
    
    % Quadratic elements (velocity)
    % Element mapping
    He=reshape(R*(P0\ze), 3,3); % zz=sum((Vq*He).*Vq,2);
    % Element Jacobian Matrix
    Je=2*He(:,2:3);
    % Element Jacobian Determinant
    jac=imag(conj(Vq*Je(:,1)).*(Vq*Je(:,2)));
    % Element Gradient
    Grad=-1i*(diag((Vq*Je(:,2))./jac)*Dx_Phi2 + ...
             -diag((Vq*Je(:,1))./jac)*Dy_Phi2 );
    % Element mass and stiffness    
    Me=Phi2'*diag(wq.*jac)*Phi2;
    Ke=real(Grad'*diag(wq.*jac)*Grad);
    [ni,nj]=ndgrid(nodes(1:m2));
    ei2(:,e)=ni(:);
    ej2(:,e)=nj(:);
    M2(:,e)=Me(:);
    K2(:,e)=Ke(:);
    

    
    % Gradient force (coupling velocity + pressure)
    Sz=Phi1'*diag(wq.*jac)*Grad;    
    [ni,nj]=ndgrid(big2small(nodes(1:m1)),nodes(1:m2));
    ei3(:,e)=ni(:);
    ej3(:,e)=nj(:);
    DX(:,e)=real(Sz(:));
    DY(:,e)=imag(Sz(:));
end

% Assembly
M1=sparse(ei1, ej1, M1, N1, N1);
K1=sparse(ei1, ej1, K1, N1, N1);
M2=sparse(ei2, ej2, M2, N2, N2);
K2=sparse(ei2, ej2, K2, N2, N2);
DX=sparse(ei3, ej3, DX, N1, N2);
DY=sparse(ei3, ej3, DY, N1, N2);

nat=[];
bnd=setdiff(bnd,nat);

mask2=1:N2;
mask2(bnd)=[];

mask1=big2small(mask2);
mask1(mask1==0)=[];

R1=sparse(mask1,1:numel(mask1),1,N1,numel(mask1));
R2=sparse(mask2,1:numel(mask2),1,N2,numel(mask2));

ei4=[1:N1-1; 1:N1-1]; ei4(2,:)=N1;
ej4=[1:N1-1; 1:N1-1];
eq4=[ones(1,N1-1); -ones(1,N1-1)];

Q1=sparse(ei4,ej4,eq4,N1,N1-1);
Q2=sparse(mask2,1:numel(mask2),1,N2,numel(mask2));

A2=Q2'*K2*Q2;
U2=chol(A2,'upper');
C1=Q1'*DX*Q2;
C2=Q1'*DY*Q2;
pcalls=0;
function [b] = uzawa(x)
    b=U2\(U2'\[C1'*x, C2'*x]);
    b=-C1*b(:,1)-C2*b(:,2);
    pcalls=pcalls+1;
end

% Exact solution
[uex,vex,pex]=stokes_exact(real(vert),imag(vert),2,1,4);

% Poisson right hand side
u0 = zeros(size(uex));
v0 = zeros(size(vex));
u0(bnd)=uex(bnd);
v0(bnd)=vex(bnd);
fu = -K2*u0;
fv = -K2*v0;
gp =  DX*u0+DY*v0;

% Solving
% A=[ A2, sparse(zeros(size(Q2,2))), -C1';
%     sparse(zeros(size(Q2,2))), A2, -C2';
%    -C1, -C2, sparse(zeros(size(Q1,2)))];
% rhs=[Q2'*fu; Q2'*fv; Q1'*gp];
% uvp=blkdiag(Q2,Q2,Q1)*(A\rhs);
% u=u0+uvp(1:N2);
% v=v0+uvp((1+N2):(2*N2));
% p=uvp((1+2*N2:N1+2*N2));

b=U2\(U2'\(Q2'*[fu,fv]));
fp=Q1'*gp+C1*b(:,1)+C2*b(:,2);
[p,flag,relres,iter,resvec]=gmres(@uzawa,fp,restart,tol,maxit,Q1'*M1*Q1);
p=Q1*p;
u=u0+Q2*(U2\(U2'\(Q2'*(fu+DX'*p))));
v=v0+Q2*(U2\(U2'\(Q2'*(fv+DY'*p))));

% Convergence plot
figure(10);
semilogy(1:size(resvec,1),resvec);

% L2 norm of the error
erru=sqrt((u-uex)'*M2*(u-uex));
errv=sqrt((v-vex)'*M2*(v-vex));
errp=sqrt((p-pex(id1))'*M1*(p-pex(id1)));

% Interpolation to cross-section
nex=1024;
xe=linspace(-6,12,nex)';
Pe=zeros(nex,4);
Pe(:,[1,2])=[xe,xe];
Pe(:,3)=0.5;
Pe(:,4)=1.0;
Pe=reshape(Pe,[],2);
[~,~,pex]=stokes_exact(Pe(:,1),Pe(:,2),2,1,4);
pex=reshape(pex,[],2);

nq=128;
xq=linspace(-6,12,nq)';
Pq=zeros(nq,4);
Pq(:,[1,2])=[xq,xq];
Pq(:,3)=0.5;
Pq(:,4)=1.0;
Pq=reshape(Pq,[],2);
DT = delaunayTriangulation([real(vert(id1)), imag(vert(id1))]);
[ti,bc] = pointLocation(DT,Pq);
pq = dot(bc',p(DT(ti,:))');
pq = reshape(pq,[],2);

figure(4);
plot(xe,pex(:,1),'r',xq,pq(:,1),'.r', xe,pex(:,2),'b',xq,pq(:,2),'.b');
title('Pressure, $p$');
legs={'Exact ($y=0.5$)','FEM ($y=0.5$)','Exact ($y=1$)','FEM ($y=1$)'};
legend(legs,'Interpreter','latex');
set(gcf,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$x$'); xlim(xe([1,end]));

% Plotting
tri2=delaunay(real(vert),imag(vert));
kd=abs(mean(vert(tri2),2)-2i)>=1;
tri2=tri2(kd,:);

figure(1);
trisurf(tri2,real(vert),imag(vert),u,u);
colormap(jet(256)); colorbar('TickLabelInterpreter','latex');
alpha(0.8);view(2);
dx=diff(xlim()); dy=diff(ylim());
pbaspect([dx,dy,min(dx,dy)]);
set(gcf,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$x$'); ylabel('$y$'); title('$u$ velocity');

figure(2);
trisurf(tri2,real(vert),imag(vert),v,v);
colormap(jet(256)); colorbar('TickLabelInterpreter','latex');
alpha(0.8);view(2);
dx=diff(xlim()); dy=diff(ylim());
pbaspect([dx,dy,min(dx,dy)]);
set(gcf,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$x$'); ylabel('$y$'); title('$v$ velocity');

pp=zeros(size(vert));
pp(id1)=p;

figure(3);
trisurf(tri(:,1:3),real(vert),imag(vert),pp,pp);
colormap(jet(256)); colorbar('TickLabelInterpreter','latex');
alpha(0.8); view(2);
dx=diff(xlim()); dy=diff(ylim());
pbaspect([dx,dy,min(dx,dy)]);
set(gcf,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$x$'); ylabel('$y$'); title('Pressure $p$');
end


function [u,v,p]=stokes_exact(xx,yy,d,R,U)
s=sqrt(d^2-R^2);
Gamma=(d+s)./(d-s);
A=-U*d/log(Gamma);
B=2*(d+s)*U/log(Gamma);
C=2*(d-s)*U/log(Gamma);
F=U/log(Gamma);
K1=xx.^2+(s+yy).^2;
K2=xx.^2+(s-yy).^2;
u=U-F*log(K1./K2)-(2./K1).*(A+F*yy).*((s+yy)+(K1./K2).*(s-yy))+...
    -(B./K1).*((s+2*yy)-(2*yy./K1).*(s+yy).^2) + ...
    -(C./K2).*((s-2*yy)+(2*yy./K2).*(s-yy).^2);
v=2*(K2-K1)./(K1.*K2).*xx.*(A+F*yy)+...
    -(2*B./K1.^2).*xx.*yy.*(s+yy)+...
    -(2*C./K2.^2).*xx.*yy.*(s-yy);
p=-(4*B./K1.^2).*xx.*(s+yy)+...
  -(4*C./K2.^2).*xx.*(s-yy)-(16*F*s./(K1.*K2)).*xx.*yy;
p=p-mean(p);
end