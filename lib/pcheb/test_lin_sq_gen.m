%% Default parameters
N   = 10;
xnh = 40;
xel = -1:0.5:1;


%% Poisson equation with constant forcing
fbv = [];
f   = [];
c   = -1;

[dof,x,U,xx,C] = solve_lin_sq_gen(N,xel,xnh,fbv,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Modified Poisson equation with constant forcing
fbv = [];
f   = [];
c   = -1;
f0  = [];
c0  = 0;
fx  = inline('4*x.*(1 - x.^2) + x.*y.^2');
cx  = 0;
fy  = inline('1 + 0*x - y.^2');
cy  = 0;
% ALTERNATIVELY:  fy = inline('0*x - y.^2');  cy = 1;

[dof,x,U,xx,C] = solve_lin_sq_gen(N,xel,xnh,fbv,f,c,f0,c0,fx,cx,fy,cy);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Laplace equation with variable boundary-values
% This is Trefethen's Problem 36 (Chapter 13).
fbv = @fbv36;

[dof,x,U,xx,C] = solve_lin_sq_gen(N,xel,xnh,fbv);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Poisson equation with variable forcing
% This is Trefethen's Problem 16 (Chapter 7).
fbv = [];
f   = @f16;
c   = 0;

[dof,x,U,xx,C] = solve_lin_sq_gen(N,xel,xnh,fbv,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)


%% Helmholtz equation with variable forcing
% This is Trefethen's Problem 17 (Chapter 7).
k   = 9;    % pseudo-wavenumber
fbv = [];
f   = @f17; % approximates a known eigenmode of the Helmholtz eqn
c   = 0;
f0  = [];
c0  = k^2;

[dof,x,U,xx,C] = solve_lin_sq_gen(N,xel,xnh,fbv,f,c,f0,c0);

quadplot_sq(x,U,xx,C),  show_order(dof,N)


