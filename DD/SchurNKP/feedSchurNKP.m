function [S,nkp,gf,kd,gb] = feedSchurNKP(S,net,A1,B1,A2,B2,C1,C2)
% Updates SchurNKP S given domain topology net (EWNS) and NKP (A1,B1,A2,B2)
% Returns preconditioned Green's function gf and constrained basis gb = kd'
m=size(A1,1);
n=size(B2,1);
kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];

% Eigenfunctions
[E1,V1,L1]=eigenfunctions(A1, A2, C1);
[E2,V2,L2]=eigenfunctions(B2, B1, C2);

% NKP wrapper
nkp=@(uu,I,J) A1(I,:)*uu*B1(J,:)'+A2(I,:)*uu*B2(J,:)';

% Green's function
[Lx,Ly]=ndgrid(L1,L2);
LL=1./(Lx+Ly);
LL((Lx==0)|(Ly==0))=0;
function uu=greenF(F,b1,b2,b0)
    uu=zeros(m,n);
    uu(rd1,rd2)=b0;
    uu(rd1,kd2)=b1;
    uu(kd1,rd2)=b2;
    uu=E1*uu*E2';
    rhs=F-E1(:,kd1)'*(A1*uu*B1'+A2*uu*B2')*E2(:,kd2);
    uu=uu+V1(:,kd1)*(LL.*(V1(kd1,kd1)'*rhs*V2(kd2,kd2)))*V2(:,kd2)';
end
gf=@(F,b1,b2,b0) greenF(F,b1,b2,b0);
kd=@(uu) reshape(E1(:,kd1)'*uu*E2(:,kd2),[],1);
gb=@(uu) E1(:,kd1)*reshape(uu,length(kd1),length(kd2))*E2(:,kd2)';

% Schur complement
blocks=size(S,1)/m;
mask=@(x,y) sparse(x,y,1,blocks,blocks);

T1=net(1:2);
T2=net(3:4);

% Building blocks
VTA1=V1(:,kd1)'*A1*E1;
VTA2=V1(:,kd1)'*A2*E1;
VTB1=V2(:,kd2)'*B1*E2;
VTB2=V2(:,kd2)'*B2*E2;

% [E W] x [E W]
[x,y]=meshgrid(T1(T1>0), T1(T1>0) );
[I,J]=meshgrid(rd1(T1>0),rd1(T1>0));
X11=LL'*(VTA1(:,I(:)).*VTA1(:,J(:)));
X12=LL'*(VTA1(:,I(:)).*VTA2(:,J(:)));
X21=LL'*(VTA2(:,I(:)).*VTA1(:,J(:)));
X22=LL'*(VTA2(:,I(:)).*VTA2(:,J(:)));
for i=1:numel(I)
    AXX=A1(I(i),J(i))*B1+A2(I(i),J(i))*B2;
    AXX=AXX-VTB1'*diag(X11(:,i))*VTB1;
    AXX=AXX-VTB1'*diag(X12(:,i))*VTB2;
    AXX=AXX-VTB2'*diag(X21(:,i))*VTB1;
    AXX=AXX-VTB2'*diag(X22(:,i))*VTB2;
    S=S+kron(mask(x(i),y(i)), AXX);
end

% [N S] x [N S]
[x,y]=meshgrid(T2(T2>0) ,T2(T2>0) );
[I,J]=meshgrid(rd2(T2>0),rd2(T2>0));
Y11=LL*(VTB1(:,I(:)).*VTB1(:,J(:)));
Y12=LL*(VTB1(:,I(:)).*VTB2(:,J(:)));
Y21=LL*(VTB2(:,I(:)).*VTB1(:,J(:)));
Y22=LL*(VTB2(:,I(:)).*VTB2(:,J(:)));
for i=1:numel(I)
    AYY=B1(I(i),J(i))*A1+B2(I(i),J(i))*A2;
    AYY=AYY-VTA1'*diag(Y11(:,i))*VTA1;
    AYY=AYY-VTA1'*diag(Y12(:,i))*VTA2;
    AYY=AYY-VTA2'*diag(Y21(:,i))*VTA1;
    AYY=AYY-VTA2'*diag(Y22(:,i))*VTA2;
    S=S+kron(mask(x(i),y(i)), AYY);
end

% [E W] x [N S]
[x,y]=meshgrid(T1(T1>0), T2(T2>0) );
[I,J]=meshgrid(rd1(T1>0),rd2(T2>0));
for i=1:numel(I)
    AXY=B1(:,J(i))*A1(I(i),:)+B2(:,J(i))*A2(I(i),:);
    AXY=AXY-VTB1'*((VTB1(:,J(i))*VTA1(:,I(i))').*LL)*VTA1;
    AXY=AXY-VTB2'*((VTB1(:,J(i))*VTA2(:,I(i))').*LL)*VTA1;
    AXY=AXY-VTB1'*((VTB2(:,J(i))*VTA1(:,I(i))').*LL)*VTA2;
    AXY=AXY-VTB2'*((VTB2(:,J(i))*VTA2(:,I(i))').*LL)*VTA2;
    S=S+kron(mask(x(i),y(i)), AXY);
end

% [N S] x [E W]
[x,y]=meshgrid(T2(T2>0), T1(T1>0) );
[I,J]=meshgrid(rd2(T2>0),rd1(T1>0));
for i=1:numel(I)
    AYX=A1(:,J(i))*B1(I(i),:)+A2(:,J(i))*B2(I(i),:);
    AYX=AYX-VTA1'*((VTA1(:,J(i))*VTB1(:,I(i))').*LL')*VTB1;
    AYX=AYX-VTA2'*((VTA1(:,J(i))*VTB2(:,I(i))').*LL')*VTB1;
    AYX=AYX-VTA1'*((VTA2(:,J(i))*VTB1(:,I(i))').*LL')*VTB2;
    AYX=AYX-VTA2'*((VTA2(:,J(i))*VTB2(:,I(i))').*LL')*VTB2;
    S=S+kron(mask(x(i),y(i)), AYX);
end
end

function [E,V,L]=eigenfunctions(K,M,C)
% Computes stiffness SK and mass SM matrices in the constrained basis E and
% constrainded eigenfunctions V and eigenvalues L
m=size(K,1);
kd=2:m-1;
rd=[1,m];

% Constrained basis
E=eye(m);
E(rd,kd)=-C(:,kd);
E(rd,:)=C(:,rd)\E(rd,:);
SK=E'*K*E;
SM=E'*M*E;

% Eigenfunctions
V=zeros(m);
[V(kd,kd),L]=eig(SK(kd,kd), SM(kd,kd), 'vector');
V(kd,kd)=bsxfun(@rdivide, V(kd,kd), sqrt(diag(V(kd,kd)'*SM(kd,kd)*V(kd,kd)))');
V(kd,kd)=bsxfun(@rdivide, V(kd,kd), sign(V(kd(1),kd)));
V(rd,kd)=E(rd,kd)*V(kd,kd);
end