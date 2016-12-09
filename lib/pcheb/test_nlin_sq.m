%%
N    = 10;
xel  = -1:0.5:1;
xnh  = 40;
f    = inline('-exp(u.^2)');
tol  = 1e-9;

[dof,x,U,xx,C] = solve_nlin_sq(N,xel,xnh,f,tol);

quadplot_sq(x,U,xx,C),  show_order(dof,N)
