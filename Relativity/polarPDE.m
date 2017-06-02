function [] = polarPDE(m,n)
if nargin<2
    n=m;
end

[Dx,x]=chebD(2*m);
A1=diag(x.^2)*Dx*Dx+diag(x)*Dx;
[A1,Dx,x]=radial(A1,Dx,x);
L2=-[0:n/2-1 -n/2:-1].^2;
y=2*pi/n*(1:n);

C=0;
F=diag(x.^2)*(zeros(m,n));

a=1;
b=0;
E1=eye(m);
B1=a*E1(1,:)+b*Dx(1,:);
b1=4*sin(5*y);

uu=ellipticfft(A1,L2,B1,C,F,b1,1);

xx=x*cos(y);
yy=x*sin(y);

figure(1);
surf(xx(:,[end,1:end]),yy(:,[end,1:end]),uu(:,[end,1:end]));
colormap(jet(256));
shading interp; 
camlight; 
axis square;
end