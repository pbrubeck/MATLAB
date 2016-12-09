%% Default parameters
n    = 40;
h    = 0.1;
c0   = -1;
xnh0 = n;
xel4 = -1:0.5:1;
xel1 = [-1,1];
xel3 = [0, 2, 5, 10];
N4   = 8;

%% Constant forcing, 16 spectral elements (M = 4 linear elements)
% Timing the matrix inversion (see command window) ...
N   = N4;
xel = xel4;
xnh = xnh0;
s   = 0;
f   = [];
c   = c0;

[dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Interpolate to a different mesh without re-computing in full
% Using mesh spacing h = 0.1 ...
% Alternatively, we can specify the linear mesh XX in full.
xnh = h;

[C,xx] = pbary(N,xel,xnh);  quadplot_sq(x,U,xx,C)

%% Re-computing, invoking bisector symmetry to halve DOF
% Valid only if the force-function shares this symmetry pattern!
N   = N4;
xel = xel4;
xnh = xnh0;
s   = +1;
f   = [];
c   = c0;

[dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Antisymmetric forcing, invoking bisector symmetry 
N   = N4;
xel = xel4;
xnh = xnh0;
s   = -1;
f   = @fnegsym;
c   = 0;

[dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Asymmetric forcing, no bisector symmetry 
% For c=0, this corresponds to Trefethen's Problem 16 (Chapter 7).
% Notice that we can create force-functions 'on the fly' ...
N   = N4;
xel = xel4;
xnh = xnh0;
s   = 0;
f   = @f16;
c   = 0;

[dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Standard collocation method (M = 1)
% PCHEB reduces to a single element when XEL is of length M+1=2.
N   = 20;
xel = xel1;
xnh = xnh0;
s   = 0;
f   = [];
c   = c0;

[dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Using asymmetric elements (M = 3) over domain [0,10]^2 
N   = 10;
xel = xel3;
xnh = xnh0;
s   = 0;
f   = [];
c   = c0;

[dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)

%% Using a different spectral order N for each of the M linear elements 
% We can still invoke symmetry!
N   = [ 6, 10, 15 ];
xel = xel3;
xnh = xnh0;
s   = 1;
f   = [];
c   = c0;

[dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f,c);

quadplot_sq(x,U,xx,C),  show_order(dof,N)
