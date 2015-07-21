function u = PoissonPDE(f, L, n)
% Solves Poisson's equation in 2D with Dirichlet boundary conditions.
t=linspace(0, L, n);
[x,y]=meshgrid(t, t);
A=f(x(1:end), y(1:end));
A=fft2(reshape(A,n,n));

i2=[0:n/2, -n/2+1:-1].^2;
[j2,k2]=meshgrid(i2, i2);
A=A./(-j2-k2);
A(1,1)=0;

u=real(ifft2(A));
surf(x,y,u,'Edgecolor','none');
view(0, 90);
colormap(jet(128));
colorbar();
end