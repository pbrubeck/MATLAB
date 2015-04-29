function [] = sphHarmY(m, l)
% Plots the Yml spherical harmonic.
theta=linspace(0, pi, 256);
phi=linspace(0, 2*pi, 256);
c=cos(theta);
s=sin(theta);

% Evaluate Legendre associated polynomial
P=Legendre(l+1);
Pml=polyDerive(P(l+1,:), m);
A=Horner(Pml,c).*(s.^m);

figure(1);
plot(theta, A);

g=sqrt((2*l+1)*factorial(l-abs(m))/(4*pi*factorial(l+abs(m))));
rho=g*A'*exp(1i*m*phi);

figure(2);
[meshPh,meshTh]=meshgrid(phi, pi/2-theta);
[x,y,z]=sph2cart(meshPh, meshTh, real(rho).^2);
surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
shading interp
end