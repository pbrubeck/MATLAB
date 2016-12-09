%% Default parameters
N   = 10;
xnh = 60;
xel = -1:0.5:1;
rho = 0.5;
f   = [];
c   = -1;

[x,U,xx,C] = solve_lin_sq_lrefine(N,xel,rho,xnh,f,c);

quadplot_sq(x,U,xx,C),  subplot(2,2,3),  drawgrid_lrefine(xel,rho)


%% Default parameters
N   = 6;
xnh = 60;
xel = [0, 0.25, 0.6, 1];
rho = [2,4,7]/10;
f   = inline('-x.*y.*exp(10 - 50*(x.^2 + y.^2))');
c   = 0;

[x,U,xx,C] = solve_lin_sq_lrefine(N,xel,rho,xnh,f,c);

quadplot_sq(x,U,xx,C),  subplot(2,2,3),  drawgrid_lrefine(xel,rho)


