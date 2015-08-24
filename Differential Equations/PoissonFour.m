function u = PoissonFour(f, L, n)
% Solves Poisson's equation in 2D with Dirichlet periodic boundary conditions.
[xx, yy]=meshgrid(linspace(0, L, n));
[i2, j2]=meshgrid([0:n/2, -n/2+1:-1].^2);
A=fft2(f(xx, yy))./(-i2-j2); A(1,1)=0;
u=real(ifft2(A));

surf(xx,yy,u,'Edgecolor','none');
view(0, 90);
colormap(jet(128));
colorbar();
end