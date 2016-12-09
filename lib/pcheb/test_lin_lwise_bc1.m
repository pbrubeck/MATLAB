%% Default parameters
xnh = 40;  % over the full-interval [-L,L]


%% Constant forcing, 3 spectral elements (M = 1), half-width L = 1
% This corresponds to 1 linear element in the half-interval [0,L]
%   and 1 square element per quadrant
N   = 20;
L   = 1;
rho = [];  % M = 1 + length(rho)
s   = 0;
f   = [];
c   = -1;

[dof,x,U,xx,C] = solve_lin_lwise_bc1(N,L,rho,xnh,s,f,c);

quadplot_lwise(x,U,xx,C),  show_order(dof,N)

%% Symmetry invoked, 27 non-uniform elements (M = 3), half-width L = 10
N   = 8;
L   = 10;
rho = [ 0.2, 0.6 ];  % M = 1 + length(rho)
s   = +1;
f   = [];
c   = -1;

[dof,x,U,xx,C] = solve_lin_lwise_bc1(N,L,rho,xnh,s,f,c);

quadplot_lwise(x,U,xx,C),  show_order(dof,N)

%% Variable symmetry forcing, 12 non-uniform elements (M = 2), L = 1
% Symmetry pattern:  f(x,y) = f(-y,-x)
N   = 12;
L   = 1;
rho = [0, 0.4, 1]; % equivalent to:  RHO = 0.4
s   = +1;
f   = inline('-(x.^2 + y.^2)');
c   = 0;

[dof,x,U,xx,C] = solve_lin_lwise_bc1(N,L,rho,xnh,s,f,c);

quadplot_lwise(x,U,xx,C),  show_order(dof,N)

%% Check by dropping assumption of symmetry
N   = 12;
L   = 1;
rho = 0.4;
s   = 0;
f   = inline('-(x.^2 + y.^2)');
c   = 0;

[dof,x,U,xx,C] = solve_lin_lwise_bc1(N,L,rho,xnh,s,f,c);

quadplot_lwise(x,U,xx,C),  show_order(dof,N)

%% Variable antisymmetric forcing
% Symmetry pattern:  f(x,y) = -f(-y,-x)
N   = 12;
L   = 1;
rho = 0.4; 
s   = -1;
f   = inline('x.^3 - y.^3');
c   = 0;

[dof,x,U,xx,C] = solve_lin_lwise_bc1(N,L,rho,xnh,s,f,c);

quadplot_lwise(x,U,xx,C),  show_order(dof,N)

