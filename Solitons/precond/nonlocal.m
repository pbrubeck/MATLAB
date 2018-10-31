function [] = nonlocal(m,n,L)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioned Newton-Krylov method
% m Gauss-Legendre-Lobatto nodes
% n Fourier modes, must be even


%a0=24.898728258915654; a1=3.407177603769343; a2=a1; P0=0^2; sigma=6;
%a0=17.170511025351768; a1=2.602758930631574; a2=a1; P0=0^2; sigma=2;
%a0=16.154363969351561; a1=2.591730322267837; a2=a1; P0=0^2; sigma=4/3;

a0=36.275423365074928; a1=3.890147311730516; a2=a1; P0=0^2; sigma=10;
lam=-1;

% Linear Hamiltonian
L(1:2)=L;
omega=0;
VR=@(r) (omega*r).^2;
[rr,th,jac,M,H,U,hshuff,J1,J2,K]=schrodpol(m,n,L(1),0,VR);
xx=rr.*cos(th);
yy=rr.*sin(th);


VS=zeros(size(jac,1),size(jac,2));
VN=zeros(size(jac,1),size(jac,2));
VP=zeros(size(jac,1),size(jac,2));

% Ansatz
function u0=ansatz(a0,a1,a2)
om=1/a1^2;
%u0=a0*lgbeam(rr,th,0,3,om);

c=sqrt(5);
q=om*c^2;
zz=acosh((xx+1i*yy)/c);
xi=real(zz);
eta=imag(zz);
u0=a0*igbeam(xi,eta,rr,3,3,3,q,om,M);
end

% Potential
function [U]=pot(psi)
    psi2=abs(psi).^2;
    U=-K(sigma,psi2);
end

% RHS
function [r]=force(lam,psi)
    z=stiff(VN,psi);
    r=[real(z(:)); imag(z(:))];
end

% Cost function
function [E]=energy(psi)
    jpsi=J1*psi*J2';
    js=abs(jpsi).^2;
    ks=jac.*K(sigma,js);
    E=real(H(psi,psi)/2-js(:)'*ks(:)/4);
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

function [ar]=afun(r)
    r=reshape(r,m,n,2);
    js=real(VP).*(J1*r(:,:,1)*J2')+...
       imag(VP).*(J1*r(:,:,2)*J2');
    
    fr=2*J1'*(VP.*(jac.*K(sigma,js)))*J2;
    ar=zeros(m,n,2);
    ar(:,:,1)=stiff(VN,r(:,:,1))-real(fr);
    ar(:,:,2)=stiff(VN,r(:,:,2))-imag(fr);
    ar=ar(:);
end

function [pu]=pfun(u)
    uu=reshape(u,m,n,2);
    pu=zeros(m,n,2);
    pu(:,:,1)=fdm2(LL(:,:,1),V1(:,:,1),V2(:,:,1),uu(:,:,1));
    pu(:,:,2)=fdm2(LL(:,:,2),V1(:,:,2),V2(:,:,2),uu(:,:,2));
    pu=pu(:);
end

%% Newton Raphson
tol=1e-10;
maxit=3;
restart=200;
function [lam,dpsi,err,flag,relres,iter,resvec]=newton(lam,psi)

    % Set potential
    VP=J1*psi*J2';
    VN=pot(VP);
    VN=jac.*(VN-lam);
       
    for k=1:2
    VS=VN;
    % Low-Rank Approximate Jacobian
    [B,sig,A]=svds(@ashuff,[n*n,m*m],2);
    sig=sqrt(diag(sig)); 
    A=reshape(A*diag(sig),[m,m,2]);
    B=reshape(B*diag(sig),[n,n,2]);
    % Fast diagonalization
    [V1(:,:,k),L1,D1]=fdm1(A(:,:,1),A(:,:,2),1:m);
    [V2(:,:,k),L2,D2]=fdm1(B(:,:,2),B(:,:,1),1:n);
    LL(:,:,k)=L1*D2.'+D1*L2.';
    end
    
    % Krylov projection solver
    r=force(lam,psi);
    [x,flag,relres,iter,resvec]=gmres(@afun,r,restart,tol,maxit,@pfun,[],pfun(r));
    
    if(P0>0)
    mu=M(psi);
    r1=[real(mu(:)); imag(mu(:))];    
    [x1,flag1,relres1,iter1,resvec1]=gmres(@afun,r1,restart,tol,maxit,@pfun,[],pfun(r1));
    x=reshape(x,[],2)*[1;1i];
    x1=reshape(x1,[],2)*[1;1i];
    y=real(((P0-psi(:)'*mu(:))/2+mu(:)'*x(:))/(mu(:)'*x1(:)));
    lam=lam+y;
    if(flag1==3)
        flag=3;
    end
    x=x-y*x1;
    x=[real(x(:)); imag(x(:))];
    end
    
    err=abs(x'*r);
    x=reshape(x,[],2)*[1;1i];
    dpsi=reshape(x,[m,n]);
end


u0=ansatz(a0,a1,a2);
u=u0;
E=energy(u);
P=real(M(u,u));
display(E);
display(P);
ii=1:m;
jj=[1:n,1];

setlatex();
figure(1);
h1=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
xlim([-L(1),L(1)]/sqrt(2));
ylim([-L(2),L(2)]/sqrt(2));
pbaspect([L([1,2]),norm(L)/sqrt(2)]);
colormap(magma(256));
colorbar();
shading interp;
axis square;
view(2);
title(num2str(E,'$E = %f$'));
drawnow;

figure(2);
h2=surf(xx(ii,jj),yy(ii,jj),angle(u(ii,jj)));
xlim([-L(1),L(1)]/sqrt(2));
ylim([-L(2),L(2)]/sqrt(2));
pbaspect([L([1,2]),norm(L)/sqrt(2)]);
caxis manual;
caxis([-pi,pi]);
colormap(hsv(256));
colorbar();
shading interp;
axis square;
view(2);
title(num2str(E,'$E = %f$'));
drawnow;

figure(3);
h3=semilogy(1:10,ones(10,1),'--*b');
title('Residual History');
drawnow;


it=0;
itnr=40;
etol=tol;
err=1;
itgmres=0;
while ( err>etol && it<itnr )
    [lam,du,err,flag,relres,iter,resvec]=newton(lam,u);
    itgmres=itgmres+length(resvec);
    u=u-du;
    display(err);
    it=it+1;
    
    P=real(M(u,u));
    E=energy(u);
    set(h1,'ZData',abs(u(ii,jj)).^2);
    title(get(1,'CurrentAxes'),num2str(E,'$E = %f$'));
    set(h2,'ZData',angle(u(ii,jj)));
    title(get(2,'CurrentAxes'),num2str(E,'$E = %f$'));
    
    set(h3,'XData',1:length(resvec));
    set(h3,'YData',resvec);
    title(get(3,'CurrentAxes'),sprintf('Newton step %d Iterations $ = %d$',it,length(resvec)))
    drawnow;
        
    if(abs(E)>1e7)
        disp('Aborting, solution blew up.');
        figure(4);
        plot(rr(:,1),real(u0(:,1)),'--b');
        xlim([0,L(1)]/sqrt(2));
        display(E);
        return
    end
end

figure(4);
plot(J1*rr(:,1),J1*real(u(:,1)),'r',J1*rr(:,1),J1*real(u0(:,1)),'--b');
xlim([0,L(1)]/sqrt(2));
display(lam);
display(E);
display(P);
display(itgmres);


% T=2*pi;
% nframes=1024;
% pbeam(T,nframes,u,xx,yy,jac,M,H,U,J1,J2,f);
end

