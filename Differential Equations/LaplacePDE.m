function u = LaplacePDE(g, L, n)
% Solves Laplace's equation in 2D subject to a boundary condition u(x,L)=g(x).
k=1:192;
y=linspace(0,L,n);
Y=sinh(pi/L*k(:)*y);
r=csch(pi*k);
c=sineSeries(g, L, 1024);
u=dst(diag(c(k).*r)*Y, n);
image([0,L], [L,0], u','CDataMapping','scaled');
colormap(jet(256));
colorbar();
end