function [rr,th,jac,M,H,U,Ht,J1,J2] = schrodpol(m,n,L,lam,VL)
% Schrodinger equation separation in polar coordinates

% Radial SEM
% Quadrature grid excludes origin!
R=L/sqrt(2);
[Dz,z]=legD(m);
[zq,wq]=gauleg(-1,1,2*m);
J1=legC(z,zq);
D=J1*Dz;

rq=R*(zq(:)+1);
jac=R*wq(:).*rq;
K1=D'*diag(jac/R^2)*D+J1'*diag(jac.*(lam+VL(rq)))*J1;
M1=J1'*diag(jac./rq.^2)*J1;
B1=J1'*diag(jac)*J1;

% Azimutal FFT
FFT=dftmtx(n)/sqrt(n);
kt=fftshift((-n/2):(n/2-1));
K2=real(FFT'*diag((2*pi/n)*kt.^2)*FFT);
M2=(2*pi/n)*eye(n);
J2=eye(n);

% Physical Domain
r=R*(z+1);
jac=repmat((2*pi/n)*jac,1,n);
[rr,th]=ndgrid(r,(2*pi/n)*(0:n-1));

% Fast diagonalization
function [V,E]=fdm(A,B,kd)
    E=zeros(m,1);
    V=eye(m);
    [V(kd,kd),E(kd)]=eig(A(kd,kd),B(kd,kd),'vector');
    s=sqrt(diag(V'*B*V));
    V=V/diag(s);
end

% Separation of variables
lt=kt.^2;
V=zeros(m,m,n);
E=zeros(m,n);
for k=1:n
    [V(:,:,k),E(:,k)]=fdm(K1+lt(k)*M1,B1,1:m);
end

function mu=mfun(u,v)
    if(nargin==1)
        mu=B1*u*M2';
    elseif(nargin==2)
        mv=B1*v*M2';
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
    u=fft(B1*u,[],2);
    for j=1:n
       u(:,j)=V(:,:,j)'*u(:,j); 
    end
    u=exp(-1i*dt*E).*u;
    for j=1:n
       u(:,j)=V(:,:,j)*u(:,j); 
    end
    u=ifft(u,[],2);
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