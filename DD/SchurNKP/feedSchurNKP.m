function [S,nkp,gf,kd,gb] = feedSchurNKP(S,net,A1,B1,A2,B2,C1,C2)
% Updates SchurNKP S given domain topology net (EWNS) and NKP (A1,B1,A2,B2)
% Returns preconditioned Green's function gf and constrained basis gb = kd'
m=size(A1,1);
n=size(B2,1);
kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];

% NKP wrapper
nkp=@(uu,I,J) A1(I,:)*uu*B1(J,:)'+A2(I,:)*uu*B2(J,:)';

% Eigenfunctions
[E1,V1,L1,D1]=eigenfunctions(A1, A2, C1);
[E2,V2,L2,D2]=eigenfunctions(B2, B1, C2);

% Green's function
[Lx,Ly]=ndgrid(L1.*D1,L2.*D2);
LL=1./((L1.*D1)*D2'+D1*(L2.*D2)');
LL((Lx==0)|(Ly==0))=0;

if false
    OP=kron(E2(:,kd2)'*B1*E2(:,kd2), E1(:,kd1)'*A1*E1(:,kd1))+...
        kron(E2(:,kd2)'*B2*E2(:,kd2), E1(:,kd1)'*A2*E1(:,kd1));
    OP2=kron(W2(kd2,kd2),W1(kd1,kd1))*diag(1./LL(:))*kron(W2(kd2,kd2)',W1(kd1,kd1)');
    disp(norm(OP-OP2,'fro'));
    
    OP2inv=kron(V2(kd2,kd2),V1(kd1,kd1))*diag(LL(:))*kron(V2(kd2,kd2)',V1(kd1,kd1)');
    disp(norm(inv(OP)-OP2inv,'fro'));
end


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
MV1=E1'*A2*V1(:,kd1);
KV1=E1'*A1*V1(:,kd1);
MV2=E2'*B1*V2(:,kd2);
KV2=E2'*B2*V2(:,kd2);

% [E W] x [E W]
[x,y]=meshgrid(T1(T1>0), T1(T1>0) );
[I,J]=meshgrid(rd1(T1>0),rd1(T1>0));
X11=(KV1(I(:),:).*KV1(J(:),:))*LL;
X12=(KV1(I(:),:).*MV1(J(:),:))*LL;
X21=(MV1(I(:),:).*KV1(J(:),:))*LL;
X22=(MV1(I(:),:).*MV1(J(:),:))*LL;
for i=1:numel(I)
    AXX=A1(I(i),J(i))*B1+A2(I(i),J(i))*B2;
    AXX=AXX-MV2*diag(X11(i,:))*MV2';
    AXX=AXX-MV2*diag(X12(i,:))*KV2';
    AXX=AXX-KV2*diag(X21(i,:))*MV2';
    AXX=AXX-KV2*diag(X22(i,:))*KV2';
    S=S+kron(mask(x(i),y(i)), (AXX+AXX')/2);
end

% [N S] x [N S]
[x,y]=meshgrid(T2(T2>0) ,T2(T2>0) );
[I,J]=meshgrid(rd2(T2>0),rd2(T2>0));
Y11=(KV2(I(:),:).*KV2(J(:),:))*LL';
Y12=(KV2(I(:),:).*MV2(J(:),:))*LL';
Y21=(MV2(I(:),:).*KV2(J(:),:))*LL';
Y22=(MV2(I(:),:).*MV2(J(:),:))*LL';
for i=1:numel(I)
    AYY=B1(I(i),J(i))*A1+B2(I(i),J(i))*A2;
    AYY=AYY-MV1*diag(Y11(i,:))*MV1';
    AYY=AYY-MV1*diag(Y12(i,:))*KV1';
    AYY=AYY-KV1*diag(Y21(i,:))*MV1';
    AYY=AYY-KV1*diag(Y22(i,:))*KV1';
    S=S+kron(mask(x(i),y(i)), (AYY+AYY')/2);
end

% [E W] x [N S]
[x,y]=meshgrid(T1(T1>0), T2(T2>0) );
[I,J]=meshgrid(rd1(T1>0),rd2(T2>0));
for i=1:numel(I)
    AXY=B1(:,J(i))*A1(I(i),:)+B2(:,J(i))*A2(I(i),:);
    AXY=AXY-MV2*((MV2(J(i),:)'*KV1(I(i),:)).*LL)*KV1';
    AXY=AXY-KV2*((MV2(J(i),:)'*MV1(I(i),:)).*LL)*KV1';
    AXY=AXY-MV2*((KV2(J(i),:)'*KV1(I(i),:)).*LL)*MV1';
    AXY=AXY-KV2*((KV2(J(i),:)'*MV1(I(i),:)).*LL)*MV1';
    %S=S+kron(mask(x(i),y(i)), AXY);
end

% [N S] x [E W]
[x,y]=meshgrid(T2(T2>0), T1(T1>0) );
[I,J]=meshgrid(rd2(T2>0),rd1(T1>0));
for i=1:numel(I)
    AYX=A1(:,J(i))*B1(I(i),:)+A2(:,J(i))*B2(I(i),:);
    AYX=AYX-MV1*((MV1(J(i),:)'*KV2(I(i),:)).*LL)*KV2';
    AYX=AYX-KV1*((MV1(J(i),:)'*MV2(I(i),:)).*LL)*KV2';
    AYX=AYX-MV1*((KV1(J(i),:)'*KV2(I(i),:)).*LL)*MV2';
    AYX=AYX-KV1*((KV1(J(i),:)'*MV2(I(i),:)).*LL)*MV2';
    S=S+kron(mask(x(i),y(i)), AYX)+kron(mask(x(i),y(i)), AYX)';
end
end

function [E,V,L,D]=eigenfunctions(K,M,C)
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

% Right eigenvectors
V=zeros(m);
[V(kd,kd),L]=eig(SK(kd,kd), SM(kd,kd), 'vector');
% Diagonal mass matrix
% D=diag(V(kd,kd)'*SM(kd,kd)*V(kd,kd));
D=sum(V(kd,kd).*(SM(kd,kd)*V(kd,kd)), 1).';

V(rd,kd)=E(rd,kd)*V(kd,kd);
V(kd,kd)=bsxfun(@rdivide, V(kd,kd), sqrt(D).');
V(kd,kd)=bsxfun(@rdivide, V(kd,kd), sign(V(kd(1),kd)));
D=sum(V(kd,kd).*(SM(kd,kd)*V(kd,kd)), 1).';

% Left eigenvectors
W=zeros(m);
W(kd,kd)=SM(kd,kd)*V(kd,kd)*diag(D);
W(rd,kd)=E(rd,kd)*W(kd,kd);
end