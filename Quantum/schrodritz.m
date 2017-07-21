function [] = schrodritz(m,n,k,p)
% Rayleigh-Ritz method for eigenvalues

k(1:2)=k;
k1=k(1);
k2=k(2);

% Differential operators
E1=eye(m);
E2=eye(n);
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
[xx,yy]=ndgrid(x,y);
xq=linspace(-1,1,1024);
yq=linspace(-1,1,1024);

% Boundary condtions
a=[1,1;1,1];
b=0*[1,-1;1,-1];
C1=diag(a(1,:))*E1([1,m],:)+diag(b(1,:))*Dx([1,m],:);
C2=diag(a(2,:))*E2([1,n],:)+diag(b(2,:))*Dy([1,n],:);


% Equation coefficients
A11=2*ones(m,n);
A21=zeros(m,n);
A12=zeros(m,n);
A22=ones(m,n);
A0=10*real(sin((xx+1i*yy).^4));

% Self-adjoint differential operator
ell=@(uu,ux,uy) Dx*(A11.*ux+A12.*uy)+(A21.*ux+A22.*uy)*Dy'+A0.*uu;
opA=@(uu) ell(uu, Dx*uu, uu*Dy');

% Preconditioner
function b=ophat(A,B,C,D,E,x,tflag)
    X=reshape(x,floor(sqrt(numel(x))),[]);
    if strcmp(tflag,'notransp')
        b=reshape(E*diag(diag(A*X.'*D).'*C)*B,[],1);
    else
        b=reshape(D*diag(C*diag(B*X.'*E))*A,[],1);
    end
end
function op2=setop(op1,A,B,C,D,E)
    op2=@(x,tflag) op1(x,tflag)+ophat(A,B,C,D,E,x,tflag);
end

ahat=@(x,tflag) ophat(1,1,A0,1,1,x,tflag);
ahat=setop(ahat,Dx,1, A11,Dx,1);
ahat=setop(ahat,1, Dy,A12,Dx,1);
ahat=setop(ahat,Dx,1, A21,1, Dy);
ahat=setop(ahat,1, Dy,A22,1, Dy);

[U,S,V]=svds(ahat, [n*n, m*m], 2);
s=diag(S);
A1=sqrt(s(1))*reshape(V(:,1),[m,m]);
A2=sqrt(s(2))*reshape(V(:,2),[m,m]);
B1=sqrt(s(1))*reshape(U(:,1),[n,n]);
B2=sqrt(s(2))*reshape(U(:,2),[n,n]);

figure(2);
subplot(2,2,1); imagesc(A1/sqrt(s(1))); title('A1'); colormap(gray(32)); colorbar;
subplot(2,2,2); imagesc(A2/sqrt(s(2))); title('A2'); colormap(gray(32)); colorbar;
subplot(2,2,3); imagesc(B1/sqrt(s(1))); title('B1'); colormap(gray(32)); colorbar;
subplot(2,2,4); imagesc(B2/sqrt(s(2))); title('B2'); colormap(gray(32)); colorbar;

% Eigenspace approximation

% Removed and Kept Degrees of freedom
rd1=[1,m];
rd2=[1,n];
kd1=2:m-1;
kd2=2:n-1;

% Give-back matrix
G1=-C1(:,rd1)\C1(:,kd1);
G2=-C2(:,rd2)\C2(:,kd2);

% Schur complement
SA1=A1(kd1,kd1)+A1(kd1,rd1)*G1;
SA2=A2(kd1,kd1)+A2(kd1,rd1)*G1;
SB1=B1(kd2,kd2)+B1(kd2,rd2)*G2;
SB2=B2(kd2,kd2)+B2(kd2,rd2)*G2;

% Generalized eigenvectors of preconditioner
V1=zeros(m,m-2);
V2=zeros(n,n-2);
[V1(kd1,:),L1]=eig(SA1,SA2,'vector');
[V2(kd2,:),L2]=eig(SB2,SB1,'vector');
V1(rd1,:)=G1*V1(kd1,:);
V2(rd2,:)=G2*V2(kd2,:);
[~,id]=sort(abs(L1),'ascend'); L1=L1(id); V1=V1(:,id);
[~,id]=sort(abs(L2),'ascend'); L2=L2(id); V2=V2(:,id);

V1=V1(:,1:k1);
V2=V2(:,1:k2);

W1=(V1(kd1,:)'*V1(kd1,:))\V1(kd1,:)'*E1(kd1,:);
W2=(V2(kd1,:)'*V2(kd1,:))\V2(kd2,:)'*E2(kd2,:);

figure(3);
subplot(1,2,1); plot(xq,interpcheb(V1,xq));
subplot(1,2,2); plot(yq,interpcheb(V2,yq));

function [F]=dimreduction(A,B,C,D,E)
    F=colkron(D,A.')*C*colkron(E,B.').';
end

Rhat=zeros(k1*k1,k2*k2);
Rhat=Rhat+dimreduction(Dx*V1, 1*V2,A11,W1*Dx,W2*1 );
Rhat=Rhat+dimreduction( 1*V1,Dy*V2,A12,W1*Dx,W2*1 );
Rhat=Rhat+dimreduction(Dx*V1, 1*V2,A21,W1*1 ,W2*Dy);
Rhat=Rhat+dimreduction( 1*V1,Dy*V2,A22,W1*1 ,W2*Dy);
Rhat=Rhat+dimreduction( 1*V1, 1*V2,A0 ,W1*1 ,W2*1 );
Rhat=reshape(Rhat,k1,k1,k2,k2);
R=zeros(k1,k2,k1*k2);
for j=1:k2
    for i=1:k1
        R(i,j,:)=reshape(squeeze(Rhat(:,i,:,j)),[],1);
    end
end
R=reshape(R,k1*k2,k1*k2);

figure(4);
imagesc(R); colormap(gray(256)); colorbar;

[uu,L]=eigs(R,p,'sm'); L=diag(L);
[~,id]=sort(abs(L),'ascend'); L=L(id); uu=uu(:,id);
display(L*4/pi^2);

uu=V1*reshape(uu(:,p),k1,k2)*V2';

figure(1);
interpplot(xq,yq,uu);
figure(5);
interpplot(xq,yq,opA(uu)-L(p)*uu);
end

% Columnwise Kronecker product
function [C]=colkron(A,B)
    C=kron(A,ones(size(B,1),1)).*kron(ones(size(A,1),1),B);
end

function []=interpplot(xq,yq,uu)
    uuu=interpcheb(interpcheb(uu,xq,1),yq,2);
    surf(xq,yq,uuu); colormap(jet(256));
    shading interp; camlight; view(2);
end