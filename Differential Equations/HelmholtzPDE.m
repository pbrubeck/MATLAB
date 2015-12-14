function A=HelmholtzPDE(n)
% Solves the Helmholtz problem in two dimensions using Chebyshev methods
k=9;
[D, x]=chebD(n);
[xx, yy]=meshgrid(x(2:end-1));
D2=D^2+(k^2/2)*eye(n); D2=D2(2:end-1, 2:end-1);

f=exp(-10*((xx-.5).^2+(yy-1).^2));
A=lyap(D2, -f);

surf(xx, yy, A, 'EdgeColor', 'none');
view(0, 90);
colormap(jet(128));
colorbar();
end