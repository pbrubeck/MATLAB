function [E] = vnlset2(n)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioner Krylov-Newton

%% Spectral element setup
rd=[];
kd=setdiff(1:n,rd);

[Dz,z,w]=legD(n);
[zq,wq]=gauleg(-1,1,2*n);
wq=wq(:);
J=legC(z,zq);
D=J*Dz;

% Physical domain
L=10;
[xx,yy]=ndgrid(L*z,L*z);
rr=hypot(yy,xx);
tt=atan2(yy,xx);
[xxq,yyq]=ndgrid(L*zq,L*zq);
rrq=hypot(yyq,xxq);

% Geometric coefficients
jac=1*(L^2)*(wq*wq');
G11=1*(wq*wq');
G12=0*(wq*wq');
G21=0*(wq*wq');
G22=1*(wq*wq');
lam=0.5;
bess=-30*besselj(1,5*rrq).^2; 

% Potential term
V=zeros(size(jac));

function U=pot(u2)
    % U = d^2/(du du*) NE
    U=-3*(J*u2*J');
end

function f=src(u2)
    % f u = d/(du*) NE
    f=-(J*u2*J');
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
    if strcmp(tflag,'transp')
        d=1;
    else 
        d=2;
    end
    Y=zeros(n*n,1);
    X=reshape(X,[n,n]);
    Y=Y+shuff(d,J,D,G11,J,D,X);
    %Y=Y+shuff(d,J,D,G12,D,J,X);
    %Y=Y+shuff(d,D,J,G21,J,D,X);
    Y=Y+shuff(d,D,J,G22,D,J,X);
    Y=Y+shuff(d,J,J,V,J,J,X);
end


%% Fast Diagonalizarion Method
function [V,L,D]=fdm1(A,B,kd)
    V=zeros(n,n);
    L=ones(n,1);
    D=ones(n,1);
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
V1=zeros(n,n);
V2=zeros(n,n);
LL=zeros(n,n);

function [f]=mass(u)
    f=J'*(jac.*(J*u*J'))*J;
end

function [f]=stiff(V0,u)
    uu=J*u*J';
    ux=D*u*J';
    uy=J*u*D';
    f=D'*(G11.*ux+G12.*uy)*J+J'*(G21.*ux+G22.*uy)*D+J'*(V0.*uu)*J;
end

function [ru]=rfun(u,tflag)
    if strcmp(tflag,'transp')
        ru=zeros(n,n);
        ru(kd,kd)=reshape(u,length(kd),length(kd));
    else 
        ru=reshape(u(kd,kd),[],1);
    end
end

function [au]=afun(u)
    au=rfun(stiff(V,rfun(u,'transp')),'notransp');
end

function [pu]=pfun(u)
    pu=rfun(fdm2(LL,V1,V2,rfun(u,'transp')),'notransp');
end

function [r]=force(u)
    F=jac.*(lam+bess+src(u.*conj(u)));
    r=rfun(stiff(F,u),'notransp');
end

function [E]=energy(u)
    V0=jac.*(lam+bess+pot(u.*conj(u))/6);
    E=real(u(:)'*reshape(stiff(V0,u),[],1))/2;
end

%% Newton Raphson
tol=1e-12;
maxit=2;
restart=100;
function [du,flag,relres,iter,resvec]=newton(r,u,ref)
    % Set potential
    V=jac.*(lam+bess+pot(u.*conj(u)));
    
    if(ref)
    % Low-Rank Approximate Jacobian
    m=2;
    [B,sig,A]=svds(@ashuff,[n*n,n*n],m);
    sig=sqrt(diag(sig)); 
    A=reshape(A*diag(sig),[n,n,m]);
    B=reshape(B*diag(sig),[n,n,m]);
    
    % Fast diagonalization
    [V1,L1,D1]=fdm1(A(:,:,1),A(:,:,2),kd);
    [V2,L2,D2]=fdm1(B(:,:,2),B(:,:,1),kd);
    LL=L1*D2.'+D1*L2.';
    end
    
    % Krylov projection solver
    [x,flag,relres,iter,resvec]=gmres(@afun,r,restart,tol,maxit,@pfun,[],r);
    du=rfun(x,'transp');
end


spin=2;
a0=(1/2)^((spin+1)/2);
a1=(2.79903);
u0=real(a0*exp(1i*spin*tt-(rr/a1).^2).*(rr.^spin));
u=u0;

V=jac.*(lam+bess+pot(u.*conj(u)));
E=energy(u);
display(E);

setlatex();
figure(1);
h1=surf(xx,yy,abs(u).^2);
colormap(magma(256));
colorbar();
shading interp;
axis square tight;
view(2);
title(num2str(E,'$E = %f$'));

figure(2);
h2=surf(xx,yy,angle(u));
colormap(hsv(256));
colorbar();
shading interp;
axis square tight;
view(2);
title(num2str(E,'$E = %f$'));

figure(3);
h3=semilogy(1:10,1:10,'--*b');
title('Residual History');

du=ones(size(u));


it=0;
itnr=40;
etol=1e-9;
while (abs(energy(du))>etol && it<itnr ) 
    ref= true ;
    [du,flag,relres,iter,resvec]=newton(force(u),u,ref);
    u=u-du;
    it=it+1;
    
    E=energy(u);
    set(h1,'ZData',abs(u).^2);
    title(get(1,'CurrentAxes'),num2str(E,'$E = %f$'))
    set(h2,'ZData',angle(u));
    title(get(2,'CurrentAxes'),num2str(E,'$E = %f$'))
    
    set(h3,'XData',1:length(resvec));
    set(h3,'YData',resvec);
    title(get(3,'CurrentAxes'),num2str(length(resvec),'Iterations $ = %d$'))
    drawnow;
end
it

k=ceil(n/2);
figure(4);
plot(xx(:,k),(u(:,k)),'r',xx(:,k),(u0(:,k)),'--b');
xlim([0,L]);
yl=ylim();
ylim([0,yl(2)]);
end