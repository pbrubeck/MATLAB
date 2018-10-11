function [] = hankdemo(m,n,L,sig)
% Radial SEM
% Quadrature grid excludes origin!
R=L/2;
[Dz,z]=legD(m);
[zq,wq]=gauleg(-1,1,2*m);
J1=legC(z,zq);
D=J1*Dz;

rq=R*(zq(:)+1);
jac=R*wq(:).*rq;
K1=D'*diag(jac/R^2)*D;
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
[rr,th]=ndgrid(r,(2*pi/n)*(0:n-1));
xx=rr.*cos(th);
yy=rr.*sin(th);
c=2;
zz=acosh((xx+1i*yy)/c);
xi=real(zz);
eta=imag(zz);
ii=1:m;
jj=[1:n,1];

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

% Initial state
omega=1;
p=3; m1=1; m2=1;
q=omega*c^2;
u0=igbeam(xi,eta,rr,p,m1,m2,q,omega,@mfun);
u=u0;

% Apply gaussian convolution
jkrm=hankbess(n,r);
u=convgauss(sig,B1,jkrm,r,u0);
%u=hank(B1,jkrm,u0);

setlatex();
figure(1);
surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
xlim([-L,L]/sqrt(2));
ylim([-L,L]/sqrt(2));
axis square;
shading interp;
colormap(magma(256));
colorbar();
view(2);
end