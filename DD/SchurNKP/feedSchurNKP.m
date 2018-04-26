function [ischur,jschur,eschur,nkp,gf]=feedSchurNKP(GID,A1,B1,A2,B2)
% Updates SchurNKP S given domain topology net (EWNS) and NKP (A1,B1,A2,B2)
% Returns preconditioned Green's function gf
m=size(A1,1);
n=size(B2,1);
kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];

% NKP wrapper
nkp=@(uu,I,J) A1(I,:)*uu*B1(J,:)'+A2(I,:)*uu*B2(J,:)';

% Eigenfunctions
[V1,L1,D1]=eigenfunctions(A1, A2);
[V2,L2,D2]=eigenfunctions(B2, B1);

% Green's function
LL=1./((L1.*D1)*D2'+D1*(L2.*D2)');
function uu=greenF(F,uu)
    rhs=F-(A1(kd1,:)*uu*B1(kd2,:)'+A2(kd1,:)*uu*B2(kd2,:)');
    uu=uu+V1(:,kd1)*(LL.*(V1(kd1,kd1)'*rhs*V2(kd2,kd2)))*V2(:,kd2)';
end
gf=@(F,uu) greenF(F,uu);

% Schur complement
ischur=zeros(m*m,16);
jschur=zeros(m*m,16);
eschur=zeros(m*m,16);
T1=rd1(any(GID(rd1,:),2));
T2=rd2(any(GID(:,rd2),1));

% Building blocks
MV1=A2*V1(:,kd1);
KV1=A1*V1(:,kd1);
MV2=B1*V2(:,kd2);
KV2=B2*V2(:,kd2);

% [E W] x [E W]
[I,J]=meshgrid(T1,T1);
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
    [x,y]=ndgrid(GID(I(i),:), GID(J(i),:));
    
    eschur(:,i)=AXX(:);
    ischur(:,i)=x(:);
    jschur(:,i)=y(:);
end

% [N S] x [N S]
[I,J]=meshgrid(T2,T2);
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
    [x,y]=ndgrid(GID(:,I(i)), GID(:,J(i)));
    
    eschur(:,i+4)=AYY(:);
    ischur(:,i+4)=x(:);
    jschur(:,i+4)=y(:);
end

% [E W] x [N S]
[I,J]=meshgrid(T1,T2);
for i=1:numel(I)
    AXY=B1(:,J(i))*A1(I(i),:)+B2(:,J(i))*A2(I(i),:);
    AXY=AXY-MV2*((MV2(J(i),:)'*KV1(I(i),:)).*LL')*KV1';
    AXY=AXY-KV2*((MV2(J(i),:)'*MV1(I(i),:)).*LL')*KV1';
    AXY=AXY-MV2*((KV2(J(i),:)'*KV1(I(i),:)).*LL')*MV1';
    AXY=AXY-KV2*((KV2(J(i),:)'*MV1(I(i),:)).*LL')*MV1';
    
    [x,y]=ndgrid(GID(I(i),:), GID(:,J(i)));
    eschur(:,i+8)=AXY(:);
    ischur(:,i+8)=x(:);
    jschur(:,i+8)=y(:);
    
    AXY=AXY';
    [x,y]=ndgrid(GID(:,J(i)), GID(I(i),:));
    eschur(:,i+12)=AXY(:);
    ischur(:,i+12)=x(:);
    jschur(:,i+12)=y(:);
end
end

function [V,L,D]=eigenfunctions(K,M,C)
% Computes stiffness SK and mass SM matrices in the constrained basis E and
% constrainded eigenfunctions V and eigenvalues L
m=size(K,1);
kd=2:m-1;
% Eigenfunctions
V=zeros(m);
[V(kd,kd),L]=eig(K(kd,kd), M(kd,kd), 'vector');
% Diagonalized mass matrix
% D=diag(V(kd,kd)'*SM(kd,kd)*V(kd,kd));
D=sum(V(kd,kd).*(M(kd,kd)*V(kd,kd)), 1).';
V(kd,kd)=bsxfun(@rdivide, V(kd,kd), sqrt(abs(D)).');
V(kd,kd)=bsxfun(@rdivide, V(kd,kd), sign(V(kd(1),kd)));
D=sign(D);
end