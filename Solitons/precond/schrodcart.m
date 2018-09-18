function [xx,yy,jac,M,H,U,Ht,J] = schrodcart(m,L,lam,VX,VY)
% Schrodinger equation separation in cartesian coordinates

% SEM
[Dz,z]=legD(m);
[zq,wq]=gauleg(-1,1,2*m);
J=legC(z,zq);
D=J*Dz;

M1=J'*diag(L*wq)*J;
M2=J'*diag(L*wq)*J;
K1=D'*diag(wq/L)*D+J'*diag(L*wq.*(lam/2+VX(L*zq)))*J;
K2=D'*diag(wq/L)*D+J'*diag(L*wq.*(lam/2+VY(L*zq)))*J;

% Physical Domain
[xx,yy]=ndgrid(L*z,L*z);
jac=L*L*wq(:)*wq(:)';

% Fast diagonalization
function [V,E]=fdm(A,B,kd)
    E=zeros(m,1);
    V=eye(m);
    [V(kd,kd),E(kd)]=eig(A(kd,kd),B(kd,kd),'vector');
    s=sqrt(diag(V'*B*V));
    V=V/diag(s);
end

% Separation of variables
[V1,E1]=fdm(K1,M1,1:m);
[V2,E2]=fdm(K2,M2,1:m);
E=E1*ones(m,1)'+ones(m,1)*E2';

function mu=mfun(u,v)
    if(nargin==1)
        mu=M1*u*M2';
    elseif(nargin==2)
        mv=M1*v*M2';
        mu=u(:)'*mv(:);
    else
        mu=0;
    end
end

function hu=hfun(u,v)
    if(nargin==1)
        hu=K1*u*M2'+M1*u*K2';
    elseif(nargin==2)
        hv=K1*v*M2'+M1*v*K2';
        hu=u(:)'*hv(:);
    else
        hu=0;
    end
end

function u=ufun(dt,u)
% u = V*exp(1i*dt*E)*V'*M*u
    u=V1'*(M1*u*M2')*V2;
    u=exp(-1i*dt*E).*u;
    u=V1*u*V2';
end

function Y=hshuff(X,tflag)
    if strcmp(tflag,'transp')
        Y=K1(:)*(M2(:)'*X(:))+M1(:)*(K2(:)'*X(:));
    else
        Y=K2(:)*(M1(:)'*X(:))+M2(:)*(K1(:)'*X(:));
    end
end

M=@mfun;
H=@hfun;
U=@ufun;
Ht=@hshuff;
end

