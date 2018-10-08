function [] = qhopoln(n,m)
% Quantum Harmonic Oscillator
% Polar coordinates
% Natural frequencies

L=20;
omega=0.1;
VL=@(r) (omega*r).^2; 

[rr,th,jac,M,H,U,hshuff,J1,J2]=schrodpol(m,n,L,0,VL);
xx=rr.*cos(th);
yy=rr.*sin(th);
ii=1:m;
jj=[1:n,1];

skw=1.4;
 
nx=3;
ny=3;
u0 = HermitePsi([zeros(nx,1);1],xx*sqrt(omega*skw)).*...
     HermitePsi([zeros(ny,1);1],yy*sqrt(omega*skw));
 
nr=8; 
l=4;
w=skw*omega*rr.^2;
c=[zeros(1,(nr-abs(l))/2),1];
u0 = LaguerreL(c,abs(l),w).*exp(-w/2).*...
     (w.^(abs(l)/2)).*exp(1i*l*th);

 
%u0 = besselj(l,omega*rr).*exp(1i*l*th);
u0 = real(u0);
u0 = u0/sqrt(M(u0,u0));



function [E]=energy(u)
    E=real(H(u,u)/M(u,u))/2;
end

lam=2*energy(u0);

% Matrix-free solver
function [au]=afun(u)
    uu=reshape(u,m,n);
    au=H(uu)-lam*(J1'*(jac.*(J1*uu*J2'))*J2);
    au=au(:);
end

function [pu]=pfun(u)
    uu=reshape(u,m,n);
    pu=fdm2(LL,V1,V2,uu);
    pu=pu(:);
end


% Matrix-free solver
V1=zeros(m,m);
V2=zeros(n,n);
LL=zeros(m,n);

% Block Shuffled jacobian
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
        Y=Y-lam*shuff(1,J2,J1,jac,J2,J1,X);
    else
        X=reshape(X,[m,m]);
        Y=Y-lam*shuff(2,J2,J1,jac,J2,J1,X);
    end
end


% Fast Diagonalization Method
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

% Newton Raphson
tol=1e-10;
maxit=3;
restart=200;
function [du,err,flag,relres,iter,resvec]=newton(u)

    % Low-Rank Approximate Jacobian
    [B,sig,A]=svds(@ashuff,[n*n,m*m],2);
    sig=sqrt(diag(sig)); 
    A=reshape(A*diag(sig),[m,m,2]);
    B=reshape(B*diag(sig),[n,n,2]);
    % Fast diagonalization
    [V1,L1,D1]=fdm1(A(:,:,1),A(:,:,2),1:m);
    [V2,L2,D2]=fdm1(B(:,:,2),B(:,:,1),1:n);
    LL=L1*D2.'+D1*L2.';

    mu=M(u);
    r=mu(:);
    
    [x,flag,relres,iter,resvec]=gmres(@afun,r,restart,tol,maxit,@pfun,[],pfun(r));
    y=real((1+u(:)'*mu(:))/(2*mu(:)'*x(:)));
    lam2=lam+y;
    r=H(u)-lam2*mu;
    r=r(:);
    
    % Krylov projection solver
    [x,flag,relres,iter,resvec]=gmres(@afun,r,restart,tol,maxit,@pfun,[],pfun(r));
    du=reshape(x,[m,n]);
    err=abs(x'*afun(x));
    
    lam=lam2;
end

u=u0/sqrt(M(u0,u0));
E=energy(u);
display(E);

setlatex();
mytitle='$z = %f$, $E/\\omega = %f$, $P = %f$';

figure(1);
h1=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
xlim([-L,L]);
ylim([-L,L]);
colormap(magma(256));
colorbar();
shading interp;
axis square;
view(2);
title(num2str(E/omega,'$E/\\omega = %f$'));

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
title(num2str(E/omega,'$E/\\omega = %f$'));

figure(3);
h3=semilogy(1:10,1:10,'--*b');
title('Residual History');

it=0;
itnr=40;
etol=1e-12;
err=1;
itgmres=0;
while ( err>etol && it<itnr ) 
    [du,err,flag,relres,iter,resvec]=newton(u);
    u=u-du;
    it=it+1;
    itgmres=itgmres+length(resvec);
    
    E=energy(u);
    P=M(u,u);

    set(h1,'ZData',abs(u(ii,jj)).^2);
    title(get(1,'CurrentAxes'),sprintf(mytitle,it,E,P));
    set(h2,'ZData',angle(u(ii,jj)));
    title(get(2,'CurrentAxes'),sprintf(mytitle,it,E,P));
    
    set(h3,'XData',1:length(resvec));
    set(h3,'YData',resvec);
    title(get(3,'CurrentAxes'),sprintf('Newton step %d Iterations $ = %d$',it,length(resvec)))
    drawnow;
    
    
    if(abs(P)>1e5)
        figure(4);
        plot(rr(:,1),real(u0(:,1)),'--b');
        xlim([0,L]);
        display(E);
        disp('Aborting, solution blew up.');
        return
    end
end

u=u/sqrt(P);
E=energy(u);
display(E);
display(itgmres);


figure(4);
plot(J1*rr(:,1),J1*real(u(:,1)),'r',J1*rr(:,1),J1*real(u0(:,1)),'--b');
%plot(rr(:,1),real(u(:,1)),'r',rr(:,1),real(u0(:,1)),'--b');
xlim([0,L]);
display(E);


t=0;
nframes=100;
T=2*pi/omega;
dt=T/nframes;
for k=1:nframes
    u=U(dt,u);
    t=t+dt;
    set(h1,'ZData',abs(u(ii,jj)).^2);
    E=energy(u)/omega;
    P=real(M(u,u));
    title(get(1,'CurrentAxes'),sprintf(mytitle,t,E,P));
    drawnow;
end

E=energy(u);
display(E);

end