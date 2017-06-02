function [uu,res] = elliptic(A1,A2,B1,B2,C,F,b1,b2,rd1,rd2,its,tol,h)
if nargin<13
    h=-1;
end
if nargin<12
    tol=1e-10;
end
if nargin<11
    its=30;
end

m=size(A1,1);
n=size(A2,1);
kd1=1:m; kd1(rd1)=[];
kd2=1:n; kd2(rd2)=[];

% Give-back matrix
G1=-B1(:,rd1)\B1(:,kd1);
G2=-B2(:,rd2)\B2(:,kd2);

% Schur complement
S1=A1(kd1,kd1)+A1(kd1,rd1)*G1;
S2=A2(kd2,kd2)+A2(kd2,rd2)*G2;

% Poincare-Steklov operator
N1=zeros(m,length(rd1)); N1(rd1,:)=inv(B1(:,rd1));
N2=zeros(n,length(rd2)); N2(rd2,:)=inv(B2(:,rd2));
P1=eye(m)-B1'/(B1*B1')*B1;
P2=eye(n)-B2'/(B2*B2')*B2;
u0=N1*sylvester(B1*B1', B2*B2', (B1*B1')*(b1*B2')+(B1*b2)*(B2*B2'))*N2';
ub=u0+(N1*b1-u0)*P2+P1*(b2*N2'-u0);

% Eigenfunctions
V1=zeros(m,length(kd1));
V2=zeros(n,length(kd2));
[V1(kd1,:),L1]=eig(S1,'vector');
[V2(kd2,:),L2]=eig(S2,'vector');
V1(rd1,:)=G1*V1(kd1,:);
V2(rd2,:)=G2*V2(kd2,:);
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2; LL(abs(LL)<1e-9)=inf;
W1=inv(V1(kd1,:));
W2=inv(V2(kd2,:));

function rhs=eqn(uu,F)
    if isfloat(C)
        rhs=A1*uu+uu*A2'+C.*uu;
    elseif ishandle(C)
        rhs=A1*uu+uu*A2'+C(uu);
    else
        rhs=A1*uu+uu*A2';
    end
    if nargin==2
        rhs=rhs-F;
    end
    rhs=rhs(kd1, kd2);
end

function uu=green(rhs)
    uu=V1*((W1*rhs*W2')./LL)*V2';
end

% Succesive Over-Relaxation
uu=ub;
um=h.*green(eqn(uu,F));
uu=uu+um;
normb=norm(eqn(ub,F),'fro');
i=1; res=norm(eqn(uu,F),'fro')/normb;
while i<its && res>=tol
    um=um+h.*green(eqn(um));
    uu=uu+um;
    res=norm(eqn(uu,F),'fro')/normb;
    i=i+1;
end

display(i);
display(res);
end