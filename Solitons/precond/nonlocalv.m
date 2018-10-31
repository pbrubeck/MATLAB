function [] = nonlocalv(m,n,L)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioned Newton-Krylov method
% m Gauss-Legendre-Lobatto nodes
% n Fourier modes, must be even


%a0=24.898728258915654; a1=3.407177603769343; a2=a1; P0=0^2; sigma=6;
%a0=17.170511025351768; a1=2.602758930631574; a2=a1; P0=0^2; sigma=2;
%a0=16.154363969351561; a1=2.591730322267837; a2=a1; P0=0^2; sigma=4/3;

a0=33; a1=4; a2=a1; P0=0^2; sigma=10;
lam=-1;

a=[a0;a1];

% Linear Hamiltonian
L(1:2)=L;
omega=0;
VR=@(r) (omega*r).^2;
[rr,th,jac,M,H,U,hshuff,J1,J2,K]=schrodpol(m,n,L(1),0,VR);
xx=rr.*cos(th);
yy=rr.*sin(th);

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

% Cost function
function [E]=energy(psi)
    jpsi=J1*psi*J2';
    js=abs(jpsi).^2;
    ks=jac.*K(sigma,js);
    E=real((H(psi,psi)-lam*M(psi,psi))/2-js(:)'*ks(:)/4);
end

cost=@(a0,a1) energy(ansatz(a0,a1,a1));

%% Newton Raphson
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
hp=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
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

it=0;
e=1e-5;
y=ones(length(a),1);
da=ones(length(a),1);
tol=1e-12;
while( abs(y'*da)>tol && it<60 )
    vars=num2cell(a);
    [y,J]=fdiff(cost,a,e,vars);
    da=-J\y;
    a=a+da;
    it=it+1;
    
    u=ansatz(vars{:});
    set(hp,'ZData',abs(u(ii,jj)).^2);
    E=real(cost(vars{:}));
    title(num2str(E,'$E = %f$'))
    drawnow;
end

display(it);
display(a);

% T=2*pi;
% nframes=1024;
% pbeam(T,nframes,u,xx,yy,jac,M,H,U,J1,J2,f);
end


function [gradf,hessf]=fdiff(f,x,e,vars)
n=length(vars);
f1=zeros(n,1);
f2=zeros(n,1);
f11=zeros(n,n);
f12=zeros(n,n);
f21=zeros(n,n);
f22=zeros(n,n);

for i=1:n
    x0=x;
    x0(i)=x0(i)+e;
    v=num2cell(x0);
    f1(i)=f(v{:});
    x0=x;
    x0(i)=x0(i)-e;
    v=num2cell(x0);
    f2(i)=f(v{:});
end
gradf=(f1-f2)/(2*e);

for i=1:n
    for j=i:n
        x0=x;
        x0(i)=x0(i)+e;
        x0(j)=x0(j)+e;
        v=num2cell(x0);
        f11(i,j)=f(v{:});
        
        x0=x;
        x0(i)=x0(i)+e;
        x0(j)=x0(j)-e;
        v=num2cell(x0);
        f12(i,j)=f(v{:});
        
        x0=x;
        x0(i)=x0(i)-e;
        x0(j)=x0(j)+e;
        v=num2cell(x0);
        f21(i,j)=f(v{:});
        
        x0=x;
        x0(i)=x0(i)-e;
        x0(j)=x0(j)-e;
        v=num2cell(x0);
        f22(i,j)=f(v{:});
        
        f11(j,i)=f11(i,j);
        f12(j,i)=f12(i,j);
        f21(j,i)=f21(i,j);
        f22(j,i)=f22(i,j);
    end
end
hessf=(f11-f12-f21+f22)/(4*e*e);
end

