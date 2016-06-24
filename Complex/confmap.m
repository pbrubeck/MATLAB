function [z] = confmap(N)
f=@(z) 1./(z);
gv=chebGrid(N);
[u, v]=ndgrid(gv);
w=u+1i*v;
err=1;

tol=1E-10;
t=0;

dt=0.5;
z=w;
z0=z/2;
while err>tol
    dz=dt.*(f(z)-w).*(z-z0)./(f(z)-f(z0));
    z0=z;
    z=z-dz;
    dt=0.5*tanh(1./abs(dz));
    dt(abs(dz)<1E-5)=0.05;
    err=max(abs(dz(:)));
    t=t+1;
end
disp(t)

mesh(real(z), imag(z), 0*real(z))
axis square; view(2); colormap([0 0 0]);
end

