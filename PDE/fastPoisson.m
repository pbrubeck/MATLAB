function [] = fastPoisson(N)
% Solves Poisson's equation in 2D with Dirichlet boundary conditions.
% Discrete Laplacian 5-point stencil
N([1 2])=N;

x=linspace(-1, 1, N(1))';
y=linspace(-1, 1, N(2));
hx=x(2)-x(1);
hy=y(2)-y(1);

[xx, yy]=ndgrid(x,y);
lx=2/hx^2*(cos(pi/(N(1)-1)*(1:N(1)-2))-1);
ly=2/hy^2*(cos(pi/(N(2)-1)*(1:N(2)-2))-1);
[Lx, Ly]=ndgrid(lx,ly); L=Lx+Ly;

% Boundary conditions
b1=[0*y; 0.2*sin(3*pi*y)]; 
b2=[0*x, (x<0).*sin(pi*x).^4];

uu=zeros(N);
uu([1 end],:)=b1;
uu(:,[1 end])=b2;

% Poisson
F=zeros(N);
E1=zeros(N(1),2); E1([2, end-1])=1/hx^2;
E2=zeros(N(2),2); E2([2, end-1])=1/hy^2;
RHS=F-E1*b1-b2*E2';

uu(2:end-1,2:end-1)=4/prod(N-1)*dst2(dst2(RHS(2:end-1,2:end-1))./L);

figure(1);
surf(xx,yy,uu); colormap(jet(256));
camlight; shading interp; axis square;
end
