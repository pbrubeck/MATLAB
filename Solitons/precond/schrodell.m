function [xi,eta,jac,M,H,U,Ht,J1,J2] = schrodell(m,n,a,b,lam)
% Schrodinger equation separation in polar coordinates

% Semifocal distance
c=sqrt(a^2-b^2);
% xi_max=2*R
R=acosh(a/c)/2;

% Radial SEM
% Quadrature grid excludes origin!
[Dz,z]=legD(m);
[zq,wq]=gauleg(-1,1,2*m); 
J1=legC(z,zq);
D=J1*Dz;

wq=wq(:);
xiq=R*(zq(:)+1);
M1=J1'*diag(R*wq)*J1;
B1=J1'*diag(R*wq.*(c*sinh(xiq)).^2)*J1;
K1=D'*diag(wq/R)*D+lam*B1;

% Azimutal FFT
FFT=dftmtx(n)/sqrt(n);
kt=fftshift((-n/2):(n/2-1));
b2=((c^2)*pi/(2*n))*[2,0,-1,zeros(1,n-5),-1,0];
M2=(2*pi/n)*eye(n);
B2=real(FFT'*toeplitz(b2)*FFT);
K2=real(FFT'*diag((2*pi/n)*kt.^2)*FFT)+lam*B2;
J2=eye(n);

% Physical Domain
r=R*(z+1);
th=(2*pi/n)*(0:n-1);
[xi,eta]=ndgrid(r,th);
jac=c^2*(2*pi*R/n)*diag(wq)*(repmat(sinh(xiq).^2,1,n)+repmat(sin(th).^2,2*m,1));

function mu=mfun(u,v)
    if(nargin==1)
        mu=B1*u*M2'+M1*u*B2';
    elseif(nargin==2)
        mv=B1*v*M2'+M1*v*B2';
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
    % TODO propagator
    u=u;
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