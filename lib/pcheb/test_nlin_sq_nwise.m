%%
N    = 6:2:24;
xel  = [-1,0,1];
xnh  = 40;
f    = inline('-exp(u)');
tol  = 1e-9;

[dof,x,U,xx,C] = solve_nlin_sq_nwise(N,xel,xnh,f,tol);

quadplot_sq(x,U,xx,C),  show_order(dof,N)
