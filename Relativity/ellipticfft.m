function [uu,res] = ellipticfft(A1,L2,B1,C,F,b1,rd1,its,tol)
if nargin<9
    tol=1e-8;
end
if nargin<8
    its=20;
end

m=size(A1,1);
kd1=1:m; kd1(rd1)=[];

% Give-back matrix
G1=-B1(:,rd1)\B1(:,kd1);
% Schur complement
S1=A1(kd1,kd1)+A1(kd1,rd1)*G1;
% Poincare-Steklov operator
N1=zeros(m,length(rd1)); N1(rd1,:)=inv(B1(:,rd1));
ub=N1*b1;

% Eigenfunctions
V1=zeros(m,length(kd1));
[V1(kd1,:),L1]=eig(S1,'vector');
V1(rd1,:)=G1*V1(kd1,:);
[L1,LL2]=ndgrid(L1,L2);
LL=L1+LL2; LL(abs(LL)<1e-9)=inf;
W1=inv(V1(kd1,:));

A2=@(uu) ifft(bsxfun(@times, fft(uu,[],2), L2),[],2);

function rhs=eqn(uu)
    if isfloat(C)
        rhs=A1*uu+A2(uu)+C.*uu-F;
    elseif ishandle(C)
        rhs=A1*uu+A2(uu)+C(uu)-F;
    else
        rhs=A1*uu+A2(uu)-F;
    end
    rhs=rhs(kd1, :);
end

function uu=green(rhs)
    uu=ifft(V1*((W1*fft(rhs,[],2))./LL),[],2);
end

% Gauss Seidel
uu=ub;
rhs=eqn(uu);
normb=norm(eqn(ub),'fro');
i=0; res=1;
while i<its && res>=tol
    uu=uu-green(rhs);
    rhs=eqn(uu);
    res=norm(rhs,'fro')/normb;
    i=i+1;
end
end