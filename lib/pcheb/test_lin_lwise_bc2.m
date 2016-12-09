%% Default parameters
xnh = 40;  % over the full-interval [-L,L]


%% Reproducing Fig6 from Pathria Pathria & Karniadakis (1995)
N   = 10;
L   = 3;
rho = 0.4; 
s   = 0;

[x,U,xx,C] = solve_lin_lwise_bc2(N,L,rho,xnh,s);

quadplot_lwise(x,U,xx,C)

%% Invoking bisector symmetry to halve DOF
% This is OK because the BCs are symmetric:  BC(-y,-x) = BC(x,y)
% If this were not the case, the code should may run but would yield a
%   different result, corresponding to symmetrized BCs.
N   = 16;
L   = 3;
rho = []; 
s   = 1;

[x,U,xx,C] = solve_lin_lwise_bc2(N,L,rho,xnh,s);

quadplot_lwise(x,U,xx,C)
