function [] = brillTest(m,n)
n(1:length(m))=n;
err=zeros(size(m));
millis=zeros(numel(m),2);
for i=1:length(m)
    [err(i), millis(i,:)]=brillSolve(m(i),n(i));
    drawnow;
end
figure(2);
semilogy(m,err);
title('Error |\Delta|_\infty');

figure(3);
plot(m,millis);
title('Time (ms)');
end

function [err,millis] = brillSolve(m,n)
millis=[0,0];
tic;
% Simulation parameters
L=8;
r0=sqrt(2)*L;

[Dx,x]=chebD(2*m);
A1=diag(x.^2)*Dx*Dx+diag(2*x)*Dx;
[A1,Dx,x]=radial(A1,Dx,x);

[Dy,y]=chebD(n); y=y';
A2=diag(1-y.^2)*Dy*Dy-diag(2*y)*Dy;

a=[1,1;0,0];
b=[1,1;1,1];

% Imposition of boundary conditions
E1=eye(m);
E2=eye(n);
B1=a(1,1)*E1(1,:)+b(1,1)*Dx(1,:);
B2=diag(a(2,:))*E2([1,end],:)+diag(b(2,:))*Dy([1,end],:);

b1=0*y+1.000068252008242;
b2=[0*x,0*x];

% Coordinate mapping
r=r0*x;
th=acos(y);
rho=r*sin(th);
z=r*cos(th);

% Test code
C=(1/2)*exp(-rho.^2-z.^2).*(1+2*rho.^2.*(-3+rho.^2+z.^2));
F=(1/20).*(1+rho.^2+z.^2).^(-5/2).*((-6)+exp(-rho.^2-z.^2) ... 
    .*(1+rho.^2+z.^2).^2.*(1+2.*rho.^2.*((-3)+rho.^2+z.^2)).*( ...
  1+10.*(1+rho.^2+z.^2).^(1/2)));
C=diag(r.^2)*C;
F=diag(r.^2)*F;

% Solution
[green,ps,kd,sc,gb]=elliptic(A1,A2,B1,B2,1,[1,n]);

ub=ps(b1,b2);
rhs=kd(F-A1*ub-ub*A2'-C.*ub);
u0=kd(ub+green(rhs));

afun=@(uu) (sc(uu)+kd(C).*uu);
pfun=@(uu) kd(green(uu));

millis(1)=toc;
tic;

[u1,res,its]=precond(afun,pfun,rhs,u0,20,2e-16);
uu=gb(u1)+ub;
millis(2)=1000*toc;

ug=1+(1/10).*(1+rho.^2+z.^2).^(-1/2);
err=norm(kd(uu-ug),'inf');

display(its);
display(res);
display(err);

figure(1);
surf(kron([-1,1],rho), z(:,[end:-1:1,1:end]), uu(:,[end:-1:1,1:end]));

colormap(jet(256));
colorbar;
shading interp; 
%camlight; 
axis square;
xlim([-L,L]);
ylim([-L,L]);
view(2);
end