function [] = spharmY(l, m, n)
% Plots the Yml spherical harmonic.
theta=linspace(0, pi, n);
phi=linspace(0, 2*pi, 2*n);
ct=cos(theta(:));
st=sin(theta(:));

% Evaluate Legendre associated polynomial
a=0; a(l+1)=1;
P=LegendreP(a,m,ct);
g=sqrt((2*l+1)*factorial(l-abs(m))/(4*pi*factorial(l+abs(m))));
rho=g*P*exp(1i*m*phi);

r=abs(real(rho));
x=r.*(st*cos(phi));
y=r.*(st*sin(phi));
z=bsxfun(@times, ct, r);

figure(1); 
h=surf(x,y,z);

set(h, 'CData', angle(rho));
alpha(0.6);
%shading interp;
lightangle(-45,30);
h.LineStyle = 'none';
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.3;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';

caxis([-pi pi]); 
colormap(hsv(256)); 
colorbar();
axis equal;
end