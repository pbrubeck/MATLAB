function []=pdesolver(n)
[D,x]=chebD(n);
D2=D*D; D2([1, end])=D2([1, end])+1;
f=16*pi^2*sin(4*pi*x); f(1)=f(1)+0.3; f(end)=f(end)-0.1;
u=D2\f;
disp(norm(D2*u-f));
figure(1); plot(x,u);

[xx, yy]=meshgrid(x); y=x';
F=50*cos(3*pi*xx).*cos(3*pi*yy); % Define differential operator
F(1,:)=F(1,:)+exp(-10*y.^2); % Add BCs
F(n,:)=F(n,:)+exp(-10*y.^2);
F(:,1)=F(:,1)+exp(-10*x.^2);
F(:,n)=F(:,n)+exp(-10*x.^2);

uu=sylvester(D2, D2', F);
disp(norm(D2*uu+uu*D2'-F, 'fro'))
figure(2); surf(xx,yy,uu); shading interp; alpha(0.6);
end