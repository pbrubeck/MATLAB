[x, w]=GaussLegendre(0,2*pi,128);
m=1;
r=linspace(-10,10,2048)';
f=exp(1i*(m*ones(size(r))*x+(2*r)*cos(x)));
J=1i^-m/(2*pi)*(f*w');
plot(r,real(J));