function uu = PoissonCheb(f, n)
% Solves Poisson's equation in 2D with Dirichlet homogenous boundary conditions.
[D, x]=chebD(n); D2=D^2; D2=D2(2:end-1,2:end-1);
[xx, yy]=ndgrid(x,x); F=f(xx,yy);
uu(2:end-1,2:end-1)=sylvester(D2, D2', F(2:end-1,2:end-1));
surf(xx,yy,uu,'Edgecolor','none');
view(0, 90); colormap(jet(128)); colorbar();
end