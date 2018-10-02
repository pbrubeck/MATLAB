function [] = vnlsetp2(m,n,L)
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
dv=0;
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
VS=zeros(size(jac,1), size(jac,2));
VN=zeros(size(jac,1), size(jac,2), 4,4);

% Physical domain
xx=rr.*cos(th);
yy=rr.*sin(th);
ii=1:m;
jj=[1:n,1];

% Ansatz
u0=(a0.^((spin+1)/2)*exp(-(xx/a1).^2-(yy/a2).^2).*...
   ((cos(ep)*xx).^2+(sin(ep)*yy).^2).^(spin/2).*...
   (cos(del)*cos(spin*th)+1i*sin(del)*sin(spin*th)));

u0=zeros(m,n,2);
u0(:,:,1)=(a0/sqrt(2))./(1+omega*rr/2).*exp( 1i*spin*th);
u0(:,:,2)=(a0/sqrt(2))./(1+omega*rr/2).*exp(-1i*spin*th);

function F=src(psi)
    psi2=sum(abs(psi).^2,3);
    F=f1(psi2);
end

function U=pot(psi)
    uu=zeros(size(psi,1),size(psi,2),2*size(psi,3));
    uu(:,:,1:2:end)=real(psi);
    uu(:,:,2:2:end)=imag(psi);
    psi2=sum(uu.^2,3);
    df1=f1(psi2);
    df2=f2(psi2);
    U=zeros(size(J1,1),size(J2,1),4,4);
    for j=1:4
        for i=1:4
        U(:,:,i,j)=(2*uu(:,:,i).*uu(:,:,j)).*df2;
        end
        U(:,:,j,j)=U(:,:,j,j)+df1;
    end
    U(:,:,1,1)=U(:,:,1,1)+dv;
    U(:,:,2,2)=U(:,:,2,2)+dv;
    U(:,:,3,3)=U(:,:,3,3)-dv;
    U(:,:,4,4)=U(:,:,4,4)-dv;
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
V1=zeros(m,m,4);
V2=zeros(n,n,4);
LL=zeros(m,n,4);

function [f]=stiff(b,u)
    f=H(u)+J1'*(b.*(J1*u*J2'))*J2;
end

function [f]=mass(b,u)
    f=J1'*(b.*(J1*u*J2'))*J2;
end

function [au]=afun(u)
    uu=reshape(u,m,n,4);
    au=zeros(m,n,4);
    for i=1:4
    au(:,:,i)=stiff(VN(:,:,i,i),uu(:,:,i));
    for j=setdiff(1:4,i)
    au(:,:,i)=au(:,:,i)+mass(VN(:,:,i,j),uu(:,:,j));
    end
    end
    au=au(:);
end

function [pu]=pfun(u)
    uu=reshape(u,m,n,4);
    pu=zeros(m,n,4);
    for i=1:4
    pu(:,:,i)=fdm2(LL(:,:,i),V1(:,:,i),V2(:,:,i),uu(:,:,i));
    end
    pu=pu(:);
end

function [r]=force(psi)
    jpsi=zeros(size(J1,1),size(J2,1),2);
    jpsi(:,:,1)=J1*psi(:,:,1)*J2';
    jpsi(:,:,2)=J1*psi(:,:,2)*J2';
    F1=jac.*(src(jpsi)+dv);
    F2=jac.*(src(jpsi)-dv);
    z1=stiff(F1,psi(:,:,1));
    z2=stiff(F2,psi(:,:,2));
    r=[real(z1(:)); imag(z1(:)); real(z2(:)); imag(z2(:))];
end

function [E]=energy(psi)
    jpsi=zeros(size(J1,1),size(J2,1),2);
    jpsi(:,:,1)=J1*psi(:,:,1)*J2';
    jpsi(:,:,2)=J1*psi(:,:,2)*J2';
    psi2=sum(abs(jpsi).^2,3);
    hpsi=zeros(size(psi));
    hpsi(:,:,1)=stiff(-dv*jac,psi(:,:,1));
    hpsi(:,:,2)=stiff( dv*jac,psi(:,:,2));
    E=real(psi(:)'*hpsi(:)+jac(:)'*(f(psi2(:))))/2;
end

%% Newton Raphson
tol=1e-10;
maxit=3;
restart=200;
function [dpsi,err,flag,relres,iter,resvec]=newton(r,psi)
    % Set potential
    jpsi=zeros(size(J1,1),size(J2,1),2);
    jpsi(:,:,1)=J1*psi(:,:,1)*J2';
    jpsi(:,:,2)=J1*psi(:,:,2)*J2';
    VN=pot(jpsi);
    for k=1:size(VN,3)
        for l=1:size(VN,4)
        VN(:,:,k,l)=jac.*VN(:,:,k,l);
        end
    end
    
    for k=1:4
    VS=VN(:,:,k,k);
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
    [x,flag,relres,iter,resvec]=gmres(@afun,r,restart,tol,maxit,@pfun,[],pfun(r));
    err=abs(x'*afun(x));
    
    x=reshape(x,[m,n,4]);
    dpsi=zeros(m,n,2);
    dpsi(:,:,1)=x(:,:,1)+1i*x(:,:,2);
    dpsi(:,:,2)=x(:,:,3)+1i*x(:,:,4);
end

u=u0;
E=energy(u);
[S1,S2,S3,S4]=stokesparams(u(:,:,1),u(:,:,2));
stksphi=atan2(S2, S1)/2;
stkschi=atan2(S3, hypot(S1,S2))/2;
display(E);

setlatex();
figure(1);
h1=surf(xx(ii,jj),yy(ii,jj),S1(ii,jj));
xlim([-L,L]);
ylim([-L,L]);
colormap(magma(256));
colorbar();
shading interp;
axis square;
view(2);
title(num2str(E,'$E = %f$'));

figure(2);
h2=surf(xx(ii,jj),yy(ii,jj),stkschi(ii,jj));
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
    [S1,S2,S3,S4]=stokesparams(u(:,:,1),u(:,:,2));
    stksphi=atan2(S2,S1)/2;
    stkschi=atan2(S3,hypot(S1,S2))/2;
    
    set(h1,'ZData',S1(ii,jj));
    title(get(1,'CurrentAxes'),num2str(E,'$E = %f$'));
    set(h2,'ZData',stkschi(ii,jj));
    title(get(2,'CurrentAxes'),num2str(E,'$E = %f$'));
    
    set(h3,'XData',1:length(resvec));
    set(h3,'YData',resvec);
    title(get(3,'CurrentAxes'),sprintf('Newton step %d Iterations $ = %d$',it,length(resvec)))
    drawnow;
    
    if(abs(E)>1e5)
        disp('Aborting, solution blew up.');
        display(E);
        
        figure(4);
        plot(rr(:,1),real(u0(:,1,1)),'--b');
        xlim([0,L]);
        display(E);
        return
    end
end

figure(4);
plot(J1*rr(:,1),J1*real(u(:,1,1)),'r',J1*rr(:,1),J1*real(u0(:,1,1)),'--b');
%plot(rr(:,1),real(u(:,1,1)),'r',rr(:,1),real(u0(:,1,1)),'--b');
xlim([0,L]);
display(E);
display(itgmres);

% T=2*pi;
% nframes=1000;
% pbeam(T,nframes,u,xx,yy,jac,M,H,U,J1,J2,f);
end