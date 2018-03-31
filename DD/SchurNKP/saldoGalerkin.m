function [stiff,mass,A1,B1,A2,B2] = saldoGalerkin(metric,Dx,x0,xx,wx,jac,g11,g12,g22)
% Stiffness matrix from Laplacian given a metric 
% Assumes oversampled coefficients
% Also returns nearest Kronecker product approximation

m=size(Dx,1);

vol=diag(wx)*( metric.*jac)*diag(wx);
C11=diag(wx)*( metric.*g22./jac)*diag(wx);
C12=diag(wx)*(-metric.*g12./jac)*diag(wx);
C22=diag(wx)*( metric.*g11./jac)*diag(wx);

E=legC(x0, xx);
ED=E*Dx;

function vv=stiffFun(uu)
    uu=reshape(uu,[m,m]);
    ux=ED*uu*E';
    uy=E*uu*ED';
    vv=ED'*(C11.*ux+C12.*uy)*E+E'*(C12.*ux+C22.*uy)*ED;
end
stiff=@stiffFun;

mass=@(uu) reshape(E'*(vol.*(E*reshape(uu,[m,m])*E'))*E,size(uu));

function B=lowrank(A,r)
    [U0,S0,V0]=svds(A,r);
    B=U0*S0*V0';
end

% This is so crucial, and I don't know why
P11=lowrank(C11,1);
P12=C12;
P22=lowrank(C22,1);

% Block-to-row permuted operator (Van Loan, 1993)
% Here we use some nifty algebra to optimize the computation
% eg. v11=(ED1'*diag(P11*diag(E2*X*E2'))*ED1);

function b=stiff_hat(x, tflag)
    X=reshape(x,floor(sqrt(numel(x))),[]).';
    b=zeros(size(X));
    if strcmp(tflag,'transp')
        b=b+ED'*diag(P11*sum(( E*X).* E, 2))*ED;
        b=b+ E'*diag(P12*sum(( E*X).*ED, 2))*ED;
        b=b+ED'*diag(P12*sum((ED*X).* E, 2))*E;
        b=b+ E'*diag(P22*sum((ED*X).*ED, 2))*E;
    else
        b=b+ E'*diag(sum((ED*X).*ED, 2).'*P11)* E;
        b=b+ED'*diag(sum((ED*X).* E, 2).'*P12)* E;
        b=b+ E'*diag(sum(( E*X).*ED, 2).'*P12)*ED;
        b=b+ED'*diag(sum(( E*X).* E, 2).'*P22)*ED;
    end
    b=reshape(b,[],1);
end

% Nearest Kronecker Product using Lanczos SVD
[B,S,A]=svds(@stiff_hat, [m*m, m*m], 2);
s=sqrt(diag(S));
A(:,1:2)=A(:,1:2)*diag(s(1:2));
B(:,1:2)=B(:,1:2)*diag(s(1:2));
A1=reshape(A(:,1),[m,m]); A1=sign(A(1,1))*(A1+A1')/2;
B1=reshape(B(:,1),[m,m]); B1=sign(B(1,1))*(B1+B1')/2;
A2=reshape(A(:,2),[m,m]); A2=sign(A(1,1))*(A2+A2')/2;
B2=reshape(B(:,2),[m,m]); B2=sign(B(1,1))*(B2+B2')/2;
end