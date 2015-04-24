function []=spherHarm(m, l)
% Plots the Yml spherical harmonic.
theta=linspace(0, pi, 100);
phi=linspace(0, 2*pi, 100);
c=cos(theta);
s=sin(theta);

% Evaluate Legendre associated polynomial
P=Legendre(l+1);
Pml=polyDerive(P(l+1,:), m);
A=Horner(Pml, c).*(s.^m);

figure(1);
plot(theta, A);

g=sqrt((2*l+1)*factorial(l-abs(m))/(4*pi*factorial(l+abs(m))));
rho=g*A'*exp(1i*m*phi);

[aux1,aux2]=meshgrid(phi, pi/2-theta);
figure(2);
[x,y,z]=sph2cart(aux1, aux2, real(rho).^2);
surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
end