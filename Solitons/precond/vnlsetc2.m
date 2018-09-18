function [] = vnlsetc2(m,L)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioned Newton-Krylov method
% Cartesian coordinates
% m Gauss-Legendre-Lobatto nodes
% n Fourier modes, must be even

% Ansatz
%spin=0; a0=(4/2)^((spin+1)/2); a1=2;
%spin=1; a0=(0.7830/2)^((spin+1)/2); a1=2.7903;
%spin=2; a0=(0.5767/2)^((spin+1)/2); a1=3.4560;
spin=2; a0=0.218858020937638; a1=3.458266856012636;

% Nonlinear potential
VN=zeros(m,m);
s=0.005;
f=@(u2) -u2/2+0*(log(s*u2+1)-s*u2)/s^2;
f1=adiff(f,1);
f2=adiff(f,2);

% Linear Hamiltonian
lam=0.5;
VX=@(x) 0*x;
VY=@(y) 0*y;
[xx,yy,jac,M,H,U,hshuff,J]=schrodcart(m,L,lam,VX,VY);
rr=hypot(yy,xx);
th=atan2(yy,xx);

zq=gauleg(-L,L,2*m);
[xq,yq]=ndgrid(zq,zq);
rq=hypot(yq,xq);
bess=-30*(besselj(1,5*rq)).^2;

function U=pot(ju)
    u2=abs(ju).^2;
    U=bess+f(u2)+u2.*(5*f1(u2)+2*ju.*f2(u2));
end

function F=src(ju)
    u2=abs(ju).^2;
    F=bess+f(u2)+u2.*f1(u2);
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
        X=reshape(X,[m,m]);
        Y=Y+shuff(1,J,J,VN,J,J,X);
    else
        X=reshape(X,[m,m]);
        Y=Y+shuff(2,J,J,VN,J,J,X);
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
V1=zeros(m,m);
V2=zeros(m,m);
LL=zeros(m,m);

function [f]=stiff(b,u)
    uu=reshape(u,m,m);
    f=H(uu)+J'*(b.*(J*uu*J'))*J;
end

function [au]=afun(u)
    au=stiff(VN,u);
    au=au(:);
end

function [pu]=pfun(u)
    uu=reshape(u,m,m);
    pu=fdm2(LL,V1,V2,uu);
    pu=pu(:);
end

function [r]=force(u)
    ju=J*u*J';
    F=jac.*src(ju);
    r=stiff(F,u);
    r=r(:);
end

function [E]=energy(psi,u)
    ju=J*psi*J';
    Vf=jac.*(bess+f(abs(ju).^2));
    hu=stiff(Vf,u);
    E=real(u(:)'*hu(:))/2;
end

%% Newton Raphson
tol=1e-12;
maxit=2;
restart=100;
function [du,flag,relres,iter,resvec]=newton(r,u,ref)
    % Set potential
    ju=J*u*J';
    VN=jac.*pot(ju);
    
    if(ref)
    % Low-Rank Approximate Jacobian
    [B,sig,A]=svds(@ashuff,[m*m,m*m],2);
    sig=sqrt(diag(sig)); 
    A=reshape(A*diag(sig),[m,m,2]);
    B=reshape(B*diag(sig),[m,m,2]);
    
    % Fast diagonalization
    [V1,L1,D1]=fdm1(A(:,:,1),A(:,:,2),1:m);
    [V2,L2,D2]=fdm1(B(:,:,2),B(:,:,1),1:m);
    LL=L1*D2.'+D1*L2.';
    end
    
    % Krylov projection solver
    [x,flag,relres,iter,resvec]=gmres(@afun,r,restart,tol,maxit,@pfun,[],r);
    du=reshape(x,[m,m]);
end

% Ansatz
u0=real(a0*exp(1i*spin*th-(rr/a1).^2).*(rr.^spin));
u=u0;

VN=jac.*pot(J*u*J');
E=energy(u,u);
display(E);

setlatex();
figure(1);
h1=surf(xx,yy,abs(u).^2);
xlim([-L,L]);
ylim([-L,L]);
colormap(magma(256));
colorbar();
shading interp;
axis square;
view(2);
title(num2str(E,'$E = %f$'));

figure(2);
h2=surf(xx,yy,angle(u));
xlim([-L,L]);
ylim([-L,L]);
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
etol=1e-9;
du=ones(size(u));
while (abs(energy(u,du))>etol && it<itnr ) 
    ref= true ;
    [du,flag,relres,iter,resvec]=newton(force(u),u,ref);
    u=u-du;
    it=it+1;
    
    E=energy(u,u)/(2*pi);
    set(h1,'ZData',abs(u).^2);
    title(get(1,'CurrentAxes'),num2str(E,'$E = %f$'))
    set(h2,'ZData',angle(u));
    title(get(2,'CurrentAxes'),num2str(E,'$E = %f$'))
    
    set(h3,'XData',1:length(resvec));
    set(h3,'YData',resvec);
    title(get(3,'CurrentAxes'),sprintf('Newton step %d Iterations $ = %d$',it,length(resvec)))
    drawnow;
end


figure(4);
k=ceil(m/2);
plot(xx(:,k),abs(u(:,k)),'r',xx(:,k),abs(u0(:,k)),'--b');
xlim([0,L]);
yl=ylim();
ylim([0,yl(2)]);
display(E);

T=2*pi;
nframes=1000;
pbeam(T,nframes,u,xx,yy,jac,M,H,U,J,J,f);
end