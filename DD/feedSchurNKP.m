function [S, nkp, gf, kd, gb] = feedSchurNKP(S, net, A1, B1, A2, B2, C1, C2)
% Updates SchurNKP S given domain topology net (EWNS) and NKP (A1,B1,A2,B2)
% Returns preconditioned Green's function gf and constrained basis gb = kd'
m=size(A1,1);
n=size(B2,1);
kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];

% Eigenfunctions
[K1,M1,E1,V1,L1]=eigenfunctions(A1, A2, C1);
[K2,M2,E2,V2,L2]=eigenfunctions(B2, B1, C2);

% NKP wrapper
nkp=@(uu,ix,iy) M1(ix,:)*uu*K2(iy,:)'+K1(ix,:)*uu*M2(iy,:)';

% Green's function
[Lx,Ly]=ndgrid(L1,L2);
LL=1./(Lx+Ly);
LL((Lx==0)|(Ly==0))=0;

function uu=greenF(F,b1,b2)
    uu=zeros(m,n);
    uu(rd1,kd2)=b1;
    uu(kd1,rd2)=b2;
    rhs=F-(K1*uu*M2'+M1*uu*K2');
    uu=E1*uu*E2'+V1(:,kd1)*(LL.*(V1(:,kd1)'*rhs*V2(:,kd2)))*V2(:,kd2)';
end
gf=@(F,b1,b2) greenF(F,b1,b2);
kd=@(uu) reshape(E1(:,kd1)'*uu*E2(:,kd2),[],1);
gb=@(uu) E1(:,kd1)*reshape(uu,length(kd1),length(kd2))*E2(:,kd2)';


blocks=size(S,1)/length(kd1);
mask=@(x,y) sparse(x,y,1,blocks,blocks);

T1=net(1:2);
T2=net(3:4);

% Building blocks
Q1=M1(kd1,kd1)*V1(kd1,kd1);
Q2=M2(kd2,kd2)*V2(kd2,kd2);
VTML1=V1(kd1,kd1)'*M1(rd1,kd1)';
VTKL1=V1(kd1,kd1)'*K1(rd1,kd1)';
VTMR1=V1(kd1,kd1)'*M1(kd1,rd1) ;
VTKR1=V1(kd1,kd1)'*K1(kd1,rd1) ;
VTML2=V2(kd2,kd2)'*M2(rd2,kd2)';
VTKL2=V2(kd2,kd2)'*K2(rd2,kd2)';
VTMR2=V2(kd2,kd2)'*M2(kd2,rd2) ;
VTKR2=V2(kd2,kd2)'*K2(kd2,rd2) ;

% S11
[x,y]=meshgrid(T1(T1>0),T1(T1>0));
[I,J]=meshgrid(find(T1>0),find(T1>0));
MM=VTML1(:,I(:)).*VTMR1(:,J(:));
MK=VTML1(:,I(:)).*VTKR1(:,J(:));
KM=VTKL1(:,I(:)).*VTMR1(:,J(:));
KK=VTKL1(:,I(:)).*VTKR1(:,J(:));
DXX=diag(L2.^2)*(LL'*MM)+diag(L2)*(LL'*(MK+KM))+(LL'*KK);
for i=1:numel(I)
    A=M1(rd1(I(i)),rd1(J(i)))*K2(kd2,kd2)+K1(rd1(I(i)),rd1(J(i)))*M2(kd2,kd2);
    S=S+kron(mask(x(i),y(i)), A-Q2*diag(DXX(:,i))*Q2');
end

% S22
[x,y]=meshgrid(T2(T2>0),T2(T2>0));
[I,J]=meshgrid(find(T2>0),find(T2>0));
MM=VTML2(:,I(:)).*VTMR2(:,J(:));
MK=VTML2(:,I(:)).*VTKR2(:,J(:));
KM=VTKL2(:,I(:)).*VTMR2(:,J(:));
KK=VTKL2(:,I(:)).*VTKR2(:,J(:));
DYY=diag(L1.^2)*(LL*MM)+diag(L1)*(LL*(MK+KM))+(LL*KK);
for i=1:numel(I)
    A=M2(rd2(I(i)),rd2(J(i)))*K1(kd1,kd1)+K2(rd2(I(i)),rd2(J(i)))*M1(kd1,kd1);
    S=S+kron(mask(x(i),y(i)), A-Q1*diag(DYY(:,i))*Q1');
end

% S12
[x,y]=meshgrid(T1(T1>0),T2(T2>0));
[I,J]=meshgrid(find(T1>0),find(T2>0));
W1=[L1.*VTML1(:,I(:)), L1.*VTKL1(:,I(:)), VTML1(:,I(:)), VTKL1(:,I(:))];
U2=[L2.*VTMR2(:,J(:)), VTMR2(:,J(:)), L2.*VTKR2(:,J(:)), VTKR2(:,J(:))];
for i=1:numel(I)
    A=M2(kd2,rd2(J(i)))*K1(rd1(I(i)),kd1)+K2(kd2,rd2(J(i)))*M1(rd1(I(i)),kd1);
    S=S+kron(mask(x(i),y(i)), A-Q2*((U2(:,i:numel(J):end)*W1(:,i:numel(I):end)').*LL')*Q1');
end

% S21
[x,y]=meshgrid(T2(T2>0),T1(T1>0));
[I,J]=meshgrid(find(T2>0),find(T1>0));
W2=[L2.*VTML2(:,I(:)), L2.*VTKL2(:,I(:)), VTML2(:,I(:)), VTKL2(:,I(:))];
U1=[L1.*VTMR1(:,J(:)), VTMR1(:,J(:)), L1.*VTKR1(:,J(:)), VTKR1(:,J(:))];
for i=1:numel(I)
    A=M1(kd1,rd1(J(i)))*K2(rd2(I(i)),kd2)+K1(kd1,rd1(J(i)))*M2(rd2(I(i)),kd2);
    S=S+kron(mask(x(i),y(i)), A-Q1*((U1(:,i:numel(J):end)*W2(:,i:numel(I):end)').*LL)*Q2');
end
end


function [SK,SM,E,V,L]=eigenfunctions(K,M,C)
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
%V(kd,kd)=bsxfun(@rdivide, V(kd,kd), sqrt(diag(V(kd,kd)'*SM(kd,kd)*V(kd,kd)))');
V(rd,kd)=E(rd,kd)*V(kd,kd);
end