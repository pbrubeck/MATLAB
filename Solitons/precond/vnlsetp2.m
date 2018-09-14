function [u] = vnlsetp2(m,n,L)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioner Krylov-Newton
% Polar coordinates
% m Gauss-Legendre-Lobatto nodes
% n Fourier modes, must be even

R=sqrt(2)*L;
lam=0.5;
bess=@(r) -30*(besselj(1,5*r)).^2;

mm=2*m;
[Dz,z,w]=legD(mm);
[zq,wq]=gauleg(-1,1,2*mm);
J=legC(z,zq);
D=J*Dz;

rq=abs(R*zq);
K=(1/R)*(D'*diag(wq.*rq)*D) + R*(J'*diag(wq.*rq.*(lam+bess(rq)))*J);
M=R*(J'*diag(wq./rq)*J);

% Symmetry BC
p1=m:-1:1;
p2=m+1:mm;
K1=K(p1,p1)+K(p1,p2)+K(p2,p1)+K(p2,p2);
M1=M(p1,p1)+M(p1,p2)+M(p2,p1)+M(p2,p2);
J1=J(mm+1:2*mm,m+1:2*m);


FFT=dftmtx(n)/(n);
kt=fftshift((-n/2):(n/2-1));
K2=real(FFT'*diag(kt.^2)*FFT);
M2=eye(n)/n;
J2=eye(n);

% Potential term
V=zeros(m,n);

% Jacobian
jac=repmat(reshape(R*wq(mm+1:2*mm).*rq(mm+1:2*mm),[],1),1,n)/n;

% Physical domain
[rr,th]=ndgrid(R*z(m+1:mm),(2*pi/n)*(0:n-1));
xx=rr.*cos(th);
yy=rr.*sin(th);

function U=pot(u2)
    % U = d^2/(du du*) NE
    U=-6*(J1*u2);
end

function f=src(u2)
    % f u = d/(du*) NE
    f=-2*(J1*u2);
end


%% Block Shuffled jacobian
function [Y]=kshuff(d,B,A,X) 
    if    (d==2)
        Y=B(:)*(A(:).'*X(:));
    elseif(d==1)
        Y=A(:)*(B(:).'*X(:));
    else
        Y=reshape(X,[],1);
    end
end

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
        Y=zeros(m*m,1);
        X=reshape(X,[n,n]);
    else
        d=2;
        Y=zeros(n*n,1);
        X=reshape(X,[m,m]);
    end
    Y=Y+kshuff(d,M2,K1,X);
    Y=Y+kshuff(d,K2,M1,X);
    Y=Y+shuff(d,J2,J1,V,J2,J1,X);
end

%% Fast Diagonalizarion Method
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
V2=zeros(n,n);
LL=zeros(m,n);

function [f]=stiff(V0,u)
    uu=reshape(u,m,n);
    f=K1*uu*M2'+M1*uu*K2'+J1'*(V0.*(J1*uu));
end

function [au]=afun(u)
    au=stiff(V,u);
    au=au(:);
end

function [pu]=pfun(u)
    uu=reshape(u,m,n);
    pu=fdm2(LL,V1,V2,uu);
    pu=pu(:);
end

function [r]=force(u)
    F=jac.*src(u.*conj(u));
    r=stiff(F,u);
    r=r(:);
end

function [E]=energy(u)
    V0=jac.*(pot(u.*conj(u))/6);
    hu=stiff(V0,u);
    E=real(u(:)'*hu(:))/2;
end

%% Newton Raphson
tol=1e-12;
maxit=2;
restart=100;
function [du,flag,relres,iter,resvec]=newton(r,u,ref)
    % Set potential
    V=jac.*pot(u.*conj(u));
    
    if(ref)
    % Low-Rank Approximate Jacobian
    [B,sig,A]=svds(@ashuff,[n*n,m*m],2);
    sig=sqrt(diag(sig)); 
    A=reshape(A*diag(sig),[m,m,2]);
    B=reshape(B*diag(sig),[n,n,2]);
    
    % Fast diagonalization
    [V1,L1,D1]=fdm1(A(:,:,1),A(:,:,2),1:m);
    [V2,L2,D2]=fdm1(B(:,:,2),B(:,:,1),1:n);
    LL=L1*D2.'+D1*L2.';
    end
    
    % Krylov projection solver
    [x,flag,relres,iter,resvec]=gmres(@afun,r,restart,tol,maxit,@pfun,[],r);
    du=reshape(x,[m,n]);
end

spin=2;
a0=(1/2)^((spin+1)/2);
a1=(2.79903);
u0=real(a0*exp(1i*spin*th-(rr/a1).^2).*(rr.^spin));
u=u0;

V=jac.*pot(u.*conj(u));
E=energy(u);
display(E);

setlatex();
figure(1);
h1=surf(xx(:,[1:end,1]),yy(:,[1:end,1]),abs(u(:,[1:end,1])).^2);
xlim([-L,L]);
ylim([-L,L]);
colormap(magma(256));
colorbar();
shading interp;
axis square;
view(2);
title(num2str(E,'$E = %f$'));

figure(2);
h2=surf(xx(:,[1:end,1]),yy(:,[1:end,1]),angle(u(:,[1:end,1])));
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
while (abs(energy(du))>etol && it<itnr ) 
    ref= true ;
    [du,flag,relres,iter,resvec]=newton(force(u),u,ref);
    u=u-du;
    it=it+1;
    
    E=energy(u);
    set(h1,'ZData',abs(u(:,[1:end,1])).^2);
    title(get(1,'CurrentAxes'),num2str(E,'$E = %f$'))
    set(h2,'ZData',angle(u(:,[1:end,1])));
    title(get(2,'CurrentAxes'),num2str(E,'$E = %f$'))
    
    set(h3,'XData',1:length(resvec));
    set(h3,'YData',resvec);
    title(get(3,'CurrentAxes'),num2str(length(resvec),'Iterations $ = %d$'))
    drawnow;
end
it

figure(4);
plot(rr(:,1),real(u(:,1)),'r',rr(:,1),real(u0(:,1)),'--b');
xlim([0,L]);
yl=ylim();
ylim([0,yl(2)]);
display(E);
end