function [gf,ps,kd,gb,dL] = elliptic(AX1,AY1,AX2,AY2,B1,B2,rd1,rd2)
m=size(AX1,1);
n=size(AY2,1);

% Kept degrees of freedom
kd1=1:m; kd1(rd1)=[];
kd2=1:n; kd2(rd2)=[];
kd=@(uu) reshape(uu(kd1,kd2),[],1);

% Give-back matrix
G1=-B1(:,rd1)\B1(:,kd1);
G2=-B2(:,rd2)\B2(:,kd2);
function uu=giveback(m,n,kd1,kd2,rd1,rd2,G1,G2,vv)
    uu=zeros(m,n);
    uu(kd1,kd2)=reshape(vv,length(kd1),length(kd2));
    uu(rd1,kd2)=G1*uu(kd1,kd2);
    uu(kd1,rd2)=uu(kd1,kd2)*G2';
    uu(rd1,rd2)=G1*uu(kd1,kd2)*G2';
end
gb=@(vv) giveback(m,n,kd1,kd2,rd1,rd2,G1,G2,vv);

% Schur complement
SX1=AX1(kd1,kd1)+AX1(kd1,rd1)*G1;
SY1=AY1(kd2,kd2)+AY1(kd2,rd2)*G2;
SX2=AX2(kd1,kd1)+AX2(kd1,rd1)*G1;
SY2=AY2(kd2,kd2)+AY2(kd2,rd2)*G2;

% Eigenfunctions
V1=zeros(m,length(kd1));
V2=zeros(n,length(kd2));
[V1(kd1,:),L1]=eig(SX1,SX2,'vector');
[V2(kd2,:),L2]=eig(SY2,SY1,'vector');
V1(rd1,:)=G1*V1(kd1,:);
V2(rd2,:)=G2*V2(kd2,:);
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2;
W1=inv(AX2(kd1,:)*V1);
W2=inv(AY1(kd2,:)*V2);

% Green's function
function uu=greenfunction(rhs,V1,V2,W1,W2)
    uu=V1*(((W1*reshape(rhs, size(V1,2), size(V2,2))*W2'))./LL)*V2';
end
gf=@(rhs) greenfunction(rhs,V1,V2,W1,W2);

% Eigenvalue perturbator
% Adds a term of the form kron(E,D)*diag(C(:))*kron(B,A)
diagop=@(A,B,C,D,E) ((A.').*D)*C*((B.').*E).';
function []=eigper(h,A,B,C,D,E,kd1,kd2,V1,V2,W1,W2)
    LL=h*LL+diagop(A*V1,B*V2,C,W1*D(kd1,:),W2*E(kd2,:));
end
dL=@(h,A,B,C,D,E) eigper(h,A,B,C,D,E,kd1,kd2,V1,V2,W1,W2); 


% Poincare-Steklov operator
N1=zeros(m,length(rd1)); N1(rd1,:)=inv(B1(:,rd1));
N2=zeros(n,length(rd2)); N2(rd2,:)=inv(B2(:,rd2));
P1=eye(m)-B1'/(B1*B1')*B1;
P2=eye(n)-B2'/(B2*B2')*B2;

u0=@(b1,b2) N1*sylvester(B1*B1', B2*B2', (B1*B1')*(b1*B2')+(B1*b2)*(B2*B2'))*N2';
ub=@(u0,b1,b2) u0+(N1*b1-u0)*P2+P1*(b2*N2'-u0);
ps=@(b1,b2) ub(u0(b1,b2),b1,b2);
end