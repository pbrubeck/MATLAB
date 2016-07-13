function uu = fastPoisson(N)
% Solves Poisson's equation in 2D with Dirichlet boundary conditions.
N([1 2])=N;
x=linspace(-1, 1, N(1))';
y=linspace(-1, 1, N(2));
hx=(x(end)-x(1))/(N(1)-1);
hy=(y(end)-y(1))/(N(2)-1);

[xx, yy]=ndgrid(x,y);
[Lx, Ly]=ndgrid(2*cos(pi/(N(1)-1)*(1:N(1)-2))-2, 2*cos(pi/(N(2)-1)*(1:N(2)-2))-2);
L=(2/hx)^2*Lx+(2/hy)^2*Ly;

% Boundary conditions
b1=[0*y; 0.2*sin(3*pi*y)];
b2=[0*x, (x<0).*sin(pi*x).^4];

% Poisson
F=zeros(N);
E1=zeros(N(1),2); E1(2,1)=1; E1(end-1,end)=1;
E2=zeros(N(2),2); E2(2,1)=1; E2(end-1,end)=1;
RHS=F-(2/hx)^2*E1*b1-b2*(2/hy)^2*E2';

uu=zeros(N);
uu([1 end],:)=b1;
uu(:,[1 end])=b2;
uu(2:end-1,2:end-1)=4/prod(N-1)*dst2(dst2(RHS(2:end-1,2:end-1))./L);

surfl(xx,yy,uu,'light');
shading interp; colormap(jet(256)); axis square;
end
