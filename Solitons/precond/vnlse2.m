function [E] = vnlse2(m,n,L)

% Ansatz
%spin=0; del=0; ep=pi/4; a0=2; a1=2; a2=a1;
%spin=1; del=pi/4; ep=pi/4; a0=1; a1=2.7903; a2=a1;
%spin=2; del=pi/4; ep=pi/4; a0=0.8291; a1=3.3941; a2=a1;
%spin=2; del=0; ep=pi/4; a0=0.5767; a1=3.4560; a2=a1;
%spin=4; del=pi/3; ep=pi/4; a0=0.421566597506070; a1=2.872534677296654; a2=a1;
%spin=2; del=0; ep=5*pi/16; a0=1.2722; a1=2.2057; a2=1.3310;

% Nonlinear potential
s=0.05;
f=@(u2) -u2/s+log(1+s*u2)/s^2;
%f=@(u2) -u2.^2/2;

% Linear Hamiltonian
lam=0.5;
VL=@(r) 0*r;
[rr,th,jac,M,H,U,hshuff,J1,J2]=schrodpol(m,n,L,lam,VL);

% Physical domain
xx=rr.*cos(th);
yy=rr.*sin(th);
ii=1:m;
jj=[1:n,1];

psi=@(a0,a1,a2) a0.^((spin+1)/2)*exp(-(xx/a1).^2-(yy/a2).^2).*...
             ((cos(ep)*xx).^2+(sin(ep)*yy).^2).^(spin/2).*...
             (cos(del)*cos(spin*th)+1i*sin(del)*sin(spin*th));

% Energy
function [E]=energy(u)
    ju=J1*u*J2';
    u2=abs(ju).^2;
    E=real(H(u,u)+jac(:)'*f(u2(:)))/2;
end

cost=@(a0,a1,a2) energy(psi(a0,a1,a2));

a=[a0;a1;a2];
vars=num2cell(a);
E=real(cost(vars{:}));
display(E);

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

it=0;
y=ones(size(a));
da=ones(size(a));
tol=1e-12;
while( abs(y'*da)>tol && it<60 )
    vars=num2cell(a);
    y=g(vars{:});
    da=-J(vars{:})\y;
    a=a+da;
    it=it+1;
    
    u=psi(vars{:});
    set(hp,'ZData',abs(u(ii,jj)).^2);
    E=real(cost(vars{:}));
    title(num2str(E,'$E = %f$'))
    drawnow;
end

display(it);
display(a);
end