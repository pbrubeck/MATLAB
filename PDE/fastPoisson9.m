function [] = fastPoisson9(N)
% Solves Poisson's equation in 2D with Dirichlet boundary conditions.
% Discrete Laplacian 9-point stencil
N([1 2])=N;
x=linspace(-1, 1, N(1))';
y=linspace(-1, 1, N(2));
hx=x(2)-x(1);
hy=y(2)-y(1);
hxy=hx*hy;

[xx, yy]=ndgrid(x,y);
lx=2/hx^2*(cos(pi/(N(1)-1)*(1:N(1)-2))-1);
ly=2/hy^2*(cos(pi/(N(2)-1)*(1:N(2)-2))-1);
L=(36-(lx*hxy-6)'*(ly*hxy-6))/(6*hxy);

% Boundary conditions
b1=[0*y; 0.2*sin(3*pi*y)];
b2=[0*x, (x<0).*sin(pi*x).^4];
uu=zeros(N);
uu([1 end],:)=b1;
uu(:,[1 end])=b2;

% Poisson
b1([1 end],[1 end])=b1([1 end],[1 end])/2;
b2([1 end],[1 end])=b2([1 end],[1 end])/2;
b1=conv2(b1, [-1, 6*hy/hx+2, -1]/(6*hxy), 'same');
b2=conv2(b2, [-1; 6*hx/hy+2; -1]/(6*hxy), 'same');

F=zeros(N);
E1=zeros(N(1),2); E1([2, end-1])=1;
E2=zeros(N(2),2); E2([2, end-1])=1;
RHS=F-E1*b1-b2*E2';

uu(2:end-1,2:end-1)=4/prod(N-1)*dst2(dst2(RHS(2:end-1,2:end-1))./L);

figure(1);
surfl(xx,yy,uu,'light'); colormap(jet(256));
shading interp; axis square;
end
