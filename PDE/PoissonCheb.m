function u = PoissonCheb(f, n)
% Solves Poisson's equation in 2D with Dirichlet homogenous boundary conditions.
[D, x]=chebD(n);
[xx, yy]=meshgrid(x(2:end-1));
D2=D^2; D2=D2(2:end-1, 2:end-1);
u=sylvester(D2, D2', f(xx, yy));

surf(xx,yy,u,'Edgecolor','none');
view(0, 90);
colormap(jet(128));
colorbar();
end