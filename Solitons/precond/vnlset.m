function [E] = vnlset(n)
% Variational Nonlinear Schrodinger Equation
% Tensor-preconditioner Krylov-Newton

%% Spectral element setup
[D,x,w]=legD(n);
hx=4;
x=hx*(x+1);
D=D/hx;
w=x.*w*hx;
K=D'*diag(w)*D;
M=diag(w);
lam=1/2;
H=K+lam*M;

% Potential term
V=zeros(size(H));

function U=pot(u2)
    % U = d^2/(du du*) NE
    U=-3*u2;
end

function f=src(u2)
    % f u = d/(du*) NE
    f=-u2;
end

%% Matrix-free solver

function [r]=force(u)
    F=diag(w.*src(u.*conj(u)));
    r=(H+F)*u;
end

function [E]=energy(u)
    E=(u'*(H+V/6)*u)/2;
end

%% Newton Raphson
tol=1e-12;

function [du]=newton(r,u)
    % Set potential
    V=diag(w.*pot(u.*conj(u)));
    du=(H+V)\r;
end

a0=sqrt(2);
a1=2;
u0=a0*exp(-(x/a1).^2);
u=u0;

V=diag(w.*pot(u.*conj(u)));
E=energy(u);
display(E);

setlatex();
figure(1);
h1=plot(x,u,'-ok');
title(num2str(E,'$E = %f$'));
drawnow;


it=0;
itnr=100;
du=ones(size(u));
while (abs(energy(du))>tol && it<itnr ) 
    du=newton(force(u),u);
    u=u-du;
    it=it+1;

    E=energy(u);
    set(h1,'YData',u);
    title(get(1,'CurrentAxes'),num2str(E,'$E = %f$'))
    drawnow;
end


figure(4);
plot(x,abs(u),'r',x,abs(u0),'--b');
end

