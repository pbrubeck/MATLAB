function u = PoissonPDE(f, L, n)
% Solves Poisson's equation in 2D with Dirichlet boundary conditions.
t=linspace(0,L,n);
[x,y]=meshgrid(t, t);
A=f(x(1:end),y(1:end));
A=fft2(reshape(A,n,n));

[j,k]=meshgrid(0:n-1, 0:n-1);
A=A./(-j.*j-k.*k);
A(1,1)=1;

u=real(ifft2(A));
image([0,L], [L,0], u','CDataMapping','scaled');
colormap(jet(64));
colorbar();
end

