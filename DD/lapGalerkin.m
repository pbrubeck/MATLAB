function [lap,mass,Mx,My,Kx,Ky] = lapGalerkin(Dx,Dy,xx,yy,wx,wy,jac,g11,g12,g22 )
% Stiffness matrix from Laplacian given a metric 
% Assumes oversampled coefficients
% Also returns nereast Kronecker product approximation

m=size(Dx,1);
n=size(Dy,1);

vol=diag(wx)*( jac)*diag(wy);
C11=diag(wx)*( g22./jac)*diag(wy);
C12=diag(wx)*(-g12./jac)*diag(wy);
C22=diag(wx)*( g11./jac)*diag(wy);

E=@(uu) interpcheb(interpcheb(uu, xx).', yy).';
E1=interpcheb(eye(m), xx);
E2=interpcheb(eye(n), yy);

function vv=stiff(uu)
    v11=Dx'*(E1'*(C11.*E(Dx*reshape(uu,[m,n])))*E2);
    v12=Dx'*(E1'*(C12.*E(reshape(uu,[m,n])))*E2)*Dy;
    v22=(E1'*(C22.*E(reshape(uu,[m,n])*Dy'))*E2)*Dy;
    vv=reshape(v11+2*v12+v22, size(uu));
end
lap=@stiff;

mass=@(uu) reshape(E1'*(vol.*E(reshape(uu,[m,n])))*E2,size(uu));

function B=lowrank(A,r)
    [U0,S0,V0]=svds(A,r);
    B=U0*S0*V0';
end

% This is sooo crucial, and I don't know why
P11=lowrank(C11,1);
P12=lowrank(C12,1);
P22=lowrank(C22,1);

function vv=laphat(x, tflag)
    X=reshape(x,floor(sqrt(numel(x))),[]);
    if strcmp(tflag,'notransp')
        v11=E2'*diag(diag(E1*(Dx*X.'*Dx')*E1').'*P11)*E2;
        v12=Dy'*E2'*diag(diag(E1*(Dx*X.')*E1').'*P12)*E2;
        v22=Dy'*E2'*diag(diag(E1*(X.')*E1').'*P22)*E2*Dy;
        vv=reshape(v11+2*v12+v22,[],1);
    else
        v22=E1'*diag(P22*diag(E2*(Dy*X*Dy')*E2'))*E1;
        v12=Dx'*E1'*diag(P12*diag(E2*(Dy*X)*E2'))*E1;
        v11=Dx'*E1'*diag(P11*diag(E2*(X)*E2'))*E1*Dx;
        vv=reshape(v11+2*v12+v22,[],1);
    end
end

[U,S,V]=svds(@laphat, [n*n, m*m], 2);
s=diag(S);
Kx=sqrt(s(1))*reshape(V(:,1),[m,m]);
Mx=sqrt(s(2))*reshape(V(:,2),[m,m]);
My=sqrt(s(1))*reshape(U(:,1),[n,n]);
Ky=sqrt(s(2))*reshape(U(:,2),[n,n]);
end