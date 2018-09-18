function [E] = vnlse2(m,n,L)

% Ansatz
%spin=0; del=0; ep=pi/4; a0=2; a1=2; a2=a1;
%spin=1; del=pi/4; ep=pi/4; a0=0.7830; a1=2.7903; a2=a1;
%spin=2; del=0; ep=pi/4; a0=0.5767; a1=3.4560; a2=a1;
%spin=4; del=pi/3; ep=pi/4; a0=0.421566597506070; a1=2.872534677296654; a2=a1;
spin=2; del=0; ep=5*pi/16; a0=1.2722; a1=2.2057; a2=1.3310;

% Nonlinear potential
s=0.05;
h=@(x) exp(-x*1E8);
g=@(x) (log(x+1)./x-1);
f=@(u2) -u2/2+0*g(s*u2)/s;

% Linear Hamiltonian
lam=0.5;
VL=@(r) -15*(besselj(0,1*r)).^2;
[rr,th,jac,M,H,U,hshuff,J1,J2]=schrodpol(m,n,L,lam,VL);

% Physical domain
xx=rr.*cos(th);
yy=rr.*sin(th);
ii=1:m;
jj=[1:n,1];

psi=@(a0,a1,a2) a0.^((spin+1)/2)*exp(-(xx/a1).^2-(yy/a2).^2).*...
             ((cos(ep)*xx).^2+(sin(ep)*yy).^2).^(spin/2).*...
             (cos(del)*cos(spin*th)+1i*sin(del)*sin(spin*th));

function [ufu]=expval(f,u)
    ju=J1*u*J2';
    u2=conj(ju).*ju;
    fu=jac.*f(u2).*ju;
    ufu=ju(:)'*fu(:);
end

% Energy
energy=@(u) (H(u,u)+expval(f,u))/2;
cost=@(a0,a1,a2) energy(psi(a0,a1,a2));

a=[a0;a1;a2];
vars=num2cell(a);
E=real(cost(vars{:}));

setlatex();
figure(1);
u=psi(vars{:});
hp=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
shading interp;
axis square;
xlim([-L,L]);
ylim([-L,L]);
colormap(magma(256));
colorbar();
view(2);
title(num2str(E,'$E = %f$'))
drawnow;

% Newton-Raphson for gradient
g=agrad(cost,length(a));
J=ahess(cost,length(a));
iter=0;
y=ones(size(a));
tol=1e-12;
while( norm(y)>tol && iter<60 )
    vars=num2cell(a);
    y=g(vars{:});
    a=a-J(vars{:})\y;
    iter=iter+1;
    
    u=psi(vars{:});
    set(hp,'ZData',abs(u(ii,jj)).^2);
    E=real(cost(vars{:}));
    title(num2str(E,'$E = %f$'))
    drawnow;
end

display(iter);
display(a);
end