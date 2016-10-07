function [] = MobiusMap( n )
r=linspace(0,1,n); r=r(:);
t=2*pi*(1:n)/n;
z=r*exp(1i*t);
w=1./z;

figure(1);
mesh(real(w),imag(w),zeros(size(w)));
shading interp; colormap([0 0 0]);
axis square; view(0,90);
end