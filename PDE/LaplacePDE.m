function u = LaplacePDE(g, L, n)
% Solves Laplace's equation in 2D subject to a boundary condition u(x,L)=g(x).
k=1:192;
y=linspace(0,pi,n);
Y=sinh(k(:)*y);
r=csch(pi*k);
c=sineSeries(g, L, 1024);
u=dst(diag(c(k).*r)*Y, n);

imagesc([0,L], [0,L], u');
set(gca, 'YDir', 'normal');
colormap(jet(256));
colorbar();
end