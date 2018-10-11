function [] = vnlset(m,n,L)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioned Newton-Krylov method
% m Gauss-Legendre-Lobatto nodes
% n Fourier modes, must be even

% Ansatz
%spin=0; del=0; ep=pi/4; a0=2; a1=2; a2=a1; P0=0;
%spin=1; del=pi/4; ep=pi/4; a0=2; a1=sqrt(8); a2=a1; P0=0;

%spin=0; del=0; ep=pi/4; a0=1.863; a1=2.529;  a2=a1; P0=12;
%spin=1; del=pi/4; ep=pi/4; a0=1.2; a1=4; a2=a1; P0=50;
%spin=2; del=pi/4; ep=pi/4; a0=1; a1=3.3941; a2=a1; P0=100;
%spin=3; del=pi/4; ep=pi/4; a0=1.65; a1=2.075; a2=a1; P0=100;

%spin=1; del=pi/4; ep=pi/4;a0=7.354225382634053; a1=0.161442173029839; a2=1; P0=110;   

spin=2; del=pi/4; ep=pi/4; a0=1; a1=3.3941; a2=a1; P0=100;

% Semifocal length
c=1;
% Nonlinear potential
s=0.05;
f=@(u2) -u2/s+log(1+s*u2)/s^2;
%f=@(u2) -u2.^2/2;
f1=adiff(f,1);
f2=adiff(f,2);

% Linear Hamiltonian
L(1:2)=L;
if (L(1)==L(2))
    omega=0;
    VR=@(r) (omega*r).^2;
    [rr,th,jac,M,H,U,hshuff,J1,J2]=schrodpol(m,n,L(1),0,VR);
    xx=rr.*cos(th);
    yy=rr.*sin(th);   
else
    [xi,eta,jac,M,H,U,hshuff,J1,J2]=schrodell(m,n,L(1),L(2),0);
    c=sqrt(L(1)^2-L(2)^2);
    xx=c*cosh(xi).*cos(eta);
    yy=c*sinh(xi).*sin(eta);
    rr=hypot(yy,xx);
    th=atan2(yy,xx);
end

VS=zeros(size(jac,1),size(jac,2));
VN=zeros(size(jac,1),size(jac,2),3);

% Ansatz
function u0=ansatz(a0,a1,a2)
u0=(a0.^((spin+1)/2)*exp(-(xx/a1).^2-(yy/a2).^2).*...
   ((cos(ep)*xx).^2+(sin(ep)*yy).^2).^(spin/2).*...
   (cos(del)*cos(spin*th)+1i*sin(del)*sin(spin*th)));
%u0=a0*igbeam(xi,eta,rr,2,2,2,a1*c^2,a1,M);
end

% Gradient
function F=src(psi)
    psi2=abs(psi).^2;
    F=f1(psi2);
end

% Hessian
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

% RHS
function [r]=force(lam,psi)
    jpsi=J1*psi*J2';
    F=jac.*(src(jpsi)-lam);
    z=stiff(F,psi);
    r=[real(z(:)); imag(z(:))];
end

% Cost function
function [E]=energy(psi)
    jpsi=J1*psi*J2';
    psi2=abs(jpsi).^2;
    E=real(H(psi,psi)+jac(:)'*f(psi2(:)))/2;
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

%% Newton Raphson
tol=1e-11;
maxit=3;
restart=200;
function [lam,dpsi,err,flag,relres,iter,resvec]=newton(lam,psi)

    % Set potential
    jpsi=J1*psi*J2';
    VN=pot(jpsi);
    VN(:,:,1)=jac.*(VN(:,:,1)-lam);
    VN(:,:,2)=jac.*(VN(:,:,2)-lam);
    VN(:,:,3)=jac.*(VN(:,:,3));
       
    for k=1:2
    VS=VN(:,:,k);
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
    x=x-y*x1;
    x=[real(x(:)); imag(x(:))];
    end
    
    err=abs(x'*r);
    x=reshape(x,[],2)*[1;1i];
    dpsi=reshape(x,[m,n]);
end

% b0=linspace(0.5,2.5,50);
% b1=linspace(0.5,2.5,50);
% ee=zeros(length(b0),length(b1));
% for i=1:length(b0)
%     for j=1:length(b1)
%         ee(i,j)=energy(ansatz(b0(i),b1(j),b1(j)));
%     end
% end
% figure(12);
% imagesc(b0,b1,ee);
% colormap(jet(256));
% colorbar();

u0=ansatz(a0,a1,a2);
u=u0;
E=energy(u);
P=real(M(u,u));
if(P0==0)
    lam=-1/2;
else
    lam=E/2;
end
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

%return;
it=0;
itnr=20;
etol=tol;
err=1;
itgmres=0;
while ( err>etol && it<itnr )
    [lam,du,err,flag,relres,iter,resvec]=newton(lam,u);
    u=u-du;
    it=it+1;
    itgmres=itgmres+length(resvec);

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
        
    if(abs(E)>1e3)
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


T=2*pi;
nframes=1024;
pbeam(T,nframes,u,xx,yy,jac,M,H,U,J1,J2,f);
end