function [] = vnlsetp(m,n,L)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioned Newton-Krylov method
% Polar coordinates
% m Gauss-Legendre-Lobatto nodes
% n Fourier modes, must be even

% Ansatz
%spin=0; del=0; ep=pi/4; a0=2; a1=2; a2=a1;
%spin=1; del=pi/4; ep=pi/4; a0=0.7830; a1=2.7903; a2=a1;
%spin=2; del=0; ep=pi/4; a0=0.5767; a1=3.4560; a2=a1;
%spin=4; del=pi/3; ep=pi/4; a0=0.421566597506070; a1=2.872534677296654; a2=a1;
%spin=2; del=0; ep=5*pi/16; a0=1.2722; a1=2.2057; a2=1.3310;

spin=2; del=0; ep=pi/4; a0=sqrt(1); a1=1; a2=1;

% Nonlinear potential
s=0.05;
f=@(u2) -u2/s+log(1+s*u2)/s^2;
%f=@(u2) -u2.^2/2;
f1=adiff(f,1);
f2=adiff(f,2);

% Linear Hamiltonian
omega=2;
lam=-omega^2;
VL=@(r) -f1((a0*besselj(spin,omega*r)).^2);
[rr,th,jac,M,H,U,hshuff,J1,J2]=schrodpol(m,n,L,lam,VL);
VS=zeros(size(jac,1),size(jac,2));
VN=zeros(size(jac,1),size(jac,2),3);

% Physical domain
xx=rr.*cos(th);
yy=rr.*sin(th);
ii=1:m;
jj=[1:n,1];

% Ansatz
u0=(a0.^((spin+1)/2)*exp(-(xx/a1).^2-(yy/a2).^2).*...
   ((cos(ep)*xx).^2+(sin(ep)*yy).^2).^(spin/2).*...
   (cos(del)*cos(spin*th)+1i*sin(del)*sin(spin*th)));

skw=1.1;
u0=a0*besselj(spin,skw*omega*rr).*exp(1i*(spin+0.01)*th);

function F=src(psi)
    psi2=abs(psi).^2;
    F=f1(psi2);
end

function U=pot(psi)
    u2=real(psi).^2;
    v2=imag(psi).^2;
    psi2=u2+v2;
    df1=f1(psi2);
    df2=f2(psi2);
    U=zeros(size(J1,1),size(J2,1),3);
    U(:,:,1)=df1+2*u2.*df2;
    U(:,:,2)=df1+2*v2.*df2;
    U(:,:,3)=2*real(psi).*imag(psi).*df2;
end


%% Block Shuffled jacobian
function [Y]=shuff(d,B2,B1,C,A2,A1,X) 
    if    (d==2)
        v=sum((B1*X).*A1,2);
        Y=reshape(B2'*diag(v'*C)*A2,[],1);
    elseif(d==1)
        v=sum((B2*X).*A2,2);
        Y=reshape(B1'*diag(C*v)*A1,[],1);
    else
        Y=reshape(X,[],1);
    end
end

function [Y]=ashuff(X,tflag)
    Y=hshuff(X,tflag);
    if strcmp(tflag,'transp')
        X=reshape(X,[n,n]);
        Y=Y+shuff(1,J2,J1,VS,J2,J1,X);
    else
        X=reshape(X,[m,m]);
        Y=Y+shuff(2,J2,J1,VS,J2,J1,X);
    end
end

%% Fast Diagonalization Method
function [V,L,D]=fdm1(A,B,kd)
    V=zeros(size(A));
    L=ones(size(A,1),1);
    D=ones(size(A,1),1);
    [V(kd,kd),L(kd)]=eig(A(kd,kd),B(kd,kd),'vector');
    D(kd)=diag(V(kd,kd)'*B(kd,kd)*V(kd,kd));
    V(kd,kd)=V(kd,kd)*diag(1./sqrt(abs(D)));
    D=sign(D);
    L=L./D;
end

function [u]=fdm2(LL,V1,V2,f)
    u=V1*((V1'*f*V2)./LL)*V2.';
end

%% Matrix-free solver
V1=zeros(m,m,2);
V2=zeros(n,n,2);
LL=zeros(m,n,2);

function [f]=stiff(b,u)
    f=H(u)+J1'*(b.*(J1*u*J2'))*J2;
end

function [f]=mass(b,u)
    f=J1'*(b.*(J1*u*J2'))*J2;
end

function [au]=afun(u)
    uu=reshape(u,m,n,2);
    au=zeros(m,n,2);
    au(:,:,1)=stiff(VN(:,:,1),uu(:,:,1))+mass(VN(:,:,3),uu(:,:,2));
    au(:,:,2)=stiff(VN(:,:,2),uu(:,:,2))+mass(VN(:,:,3),uu(:,:,1));
    au=au(:);
end

function [pu]=pfun(u)
    uu=reshape(u,m,n,2);
    pu=zeros(m,n,2);
    pu(:,:,1)=fdm2(LL(:,:,1),V1(:,:,1),V2(:,:,1),uu(:,:,1));
    pu(:,:,2)=fdm2(LL(:,:,2),V1(:,:,2),V2(:,:,2),uu(:,:,2));
    pu=pu(:);
end

function [r]=force(psi)
    jpsi=J1*psi*J2';
    F=jac.*src(jpsi);
    z=stiff(F,psi);
    r=[real(z(:)); imag(z(:))];
end

function [E]=energy(u)
    ju=J1*u*J2';
    u2=abs(ju).^2;
    Vf=jac.*f1(u2);
    hu=stiff(Vf,u);
    E=real(u(:)'*hu(:))/2;
end

%% Newton Raphson
tol=1e-10;
maxit=3;
restart=200;
function [dpsi,err,flag,relres,iter,resvec]=newton(r,psi)
    % Set potential
    jpsi=J1*psi*J2';
    VN=pot(jpsi);
    for k=1:size(VN,3)
        VN(:,:,k)=jac.*VN(:,:,k);
    end
    
    VS=VN(:,:,1);
    % Low-Rank Approximate Jacobian
    [B,sig,A]=svds(@ashuff,[n*n,m*m],2);
    sig=sqrt(diag(sig)); 
    A=reshape(A*diag(sig),[m,m,2]);
    B=reshape(B*diag(sig),[n,n,2]);
    % Fast diagonalization
    [V1(:,:,1),L1,D1]=fdm1(A(:,:,1),A(:,:,2),1:m);
    [V2(:,:,1),L2,D2]=fdm1(B(:,:,2),B(:,:,1),1:n);
    LL(:,:,1)=L1*D2.'+D1*L2.';
    
    VS=VN(:,:,2);
    % Low-Rank Approximate Jacobian
    [B,sig,A]=svds(@ashuff,[n*n,m*m],2);
    sig=sqrt(diag(sig)); 
    A=reshape(A*diag(sig),[m,m,2]);
    B=reshape(B*diag(sig),[n,n,2]);
    % Fast diagonalization
    [V1(:,:,2),L1,D1]=fdm1(A(:,:,1),A(:,:,2),1:m);
    [V2(:,:,2),L2,D2]=fdm1(B(:,:,2),B(:,:,1),1:n);
    LL(:,:,2)=L1*D2.'+D1*L2.';
    
    % Krylov projection solver
    [x,flag,relres,iter,resvec]=gmres(@afun,r,restart,tol,maxit,@pfun,[],r);
    err=abs(x'*afun(x));
    
    x=reshape(x,[m,n,2]);
    dpsi=x(:,:,1)+1i*x(:,:,2);
end

u=u0;
E=energy(u);
display(E);

setlatex();
figure(1);
h1=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
xlim([-L,L]);
ylim([-L,L]);
colormap(magma(256));
colorbar();
shading interp;
axis square;
view(2);
title(num2str(E,'$E = %f$'));

figure(2);
h2=surf(xx(ii,jj),yy(ii,jj),angle(u(ii,jj)));
xlim([-L,L]);
ylim([-L,L]);
caxis manual;
caxis([-pi,pi]);
colormap(hsv(256));
colorbar();
shading interp;
axis square;
view(2);
title(num2str(E,'$E = %f$'));

figure(3);
h3=semilogy(1:10,1:10,'--*b');
title('Residual History');

it=0;
itnr=40;
etol=1e-13;
err=1;
itgmres=0;
while ( err>etol && it<itnr )
    [du,err,flag,relres,iter,resvec]=newton(force(u),u);
    u=u-du;
    it=it+1;
    itgmres=itgmres+length(resvec);
    
    E=energy(u);
    set(h1,'ZData',abs(u(ii,jj)).^2);
    title(get(1,'CurrentAxes'),num2str(E,'$E = %f$'));
    set(h2,'ZData',angle(u(ii,jj)));
    title(get(2,'CurrentAxes'),num2str(E,'$E = %f$'));
    
    set(h3,'XData',1:length(resvec));
    set(h3,'YData',resvec);
    title(get(3,'CurrentAxes'),sprintf('Newton step %d Iterations $ = %d$',it,length(resvec)))
    drawnow;
    
    if(abs(E)>1e5)
        disp('Aborting, solution blew up.');
        display(E);
        
        figure(4);
        plot(rr(:,1),real(u0(:,1)),'--b');
        xlim([0,L]);
        display(E);
        return
    end
end

figure(4);
plot(J1*rr(:,1),J1*real(u(:,1)),'r',J1*rr(:,1),J1*real(u0(:,1)),'--b');
%plot(rr(:,1),real(u(:,1)),'r',rr(:,1),real(u0(:,1)),'--b');
xlim([0,L]);
display(E);
display(itgmres);

T=2*pi;
nframes=1000;
pbeam(T,nframes,u,xx,yy,jac,M,H,U,J1,J2,f1);
end