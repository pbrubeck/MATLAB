function [E] = vnlset1(m)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioner Krylov-Newton

%% Spectral element setup
L=10;
lam=0.5;
bess=@(r) -0*(besselj(1,1*r)).^2;
R=L/(2);

mm=m;
[Dz,z]=legD(mm);
[zq,wq]=gauleg(-1,1,2*mm);
J=legC(z,zq);
D=J*Dz;

z=z+1;
zq=zq+1;
rq=abs(R*zq);
H0=(1/R)*(D'*diag(wq.*rq)*D) + R*(J'*diag(wq.*rq.*(lam+bess(rq)))*J);

% Symmetry BC
p1=m:-1:1;
p2=m+1:mm;
p3=mm+1:2*mm;
%H0=H0(p1,p1)+H0(p1,p2)+H0(p2,p1)+H0(p2,p2);
%J=J(p3,p2);

jac=reshape(R*wq.*rq,[],1);
%jac=jac(p3);
x=abs(R*z);

% Potential term
s=0.005;
f=@(u2) -u2/2+0*(log(s*u2+1)-s*u2)/s^2;
f1=adiff(f,1);
f2=adiff(f,2);

function [E]=energy(u,psi)
    u2=abs(J*u).^2;
    V=J'*diag(jac.*f(u2))*J;
    E=(psi'*(H0+V)*psi)/2;
end

function [H]=nhess(u)
    ju=J*u;
    u2=abs(ju).^2;
    H = H0 + J'*diag(jac.*(f(u2)+u2.*(5*f1(u2)+2*ju.*f2(u2))))*J;
end

function [F]=ngrad(u)
    u2=abs(J*u).^2;
    F = H0 + J'*diag(jac.*(f(u2)+u2.*f1(u2)))*J;
end

%% Newton Raphson
tol=1e-12;

function [du]=newton(u)
    du=nhess(u)\(ngrad(u)*u);
end

a0=sqrt(2);
a1=2;
u0=a0*exp(-(x/a1).^2);
u=u0;

E=energy(u,u);
display(E);

setlatex();
figure(1);
h1=plot(x,u,'-ok');
title(num2str(E,'$E = %f$'));
drawnow;

it=0;
itnr=100;
du=ones(size(u));
while (abs(energy(u,du))>tol && it<itnr ) 
    du=newton(u);
    u=u-du;
    it=it+1;

    E=energy(u,u);
    set(h1,'YData',u);
    title(get(1,'CurrentAxes'),num2str(E,'$E = %f$'))
    drawnow;
end

figure(4);
plot(x,abs(u),'r',x,abs(u0),'--b');
end

