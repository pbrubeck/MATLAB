%% Square domain, constant forcing
N  = 8;
el = -1:0.5:1;
nh = 40;
f  = [];
c  = -1;

[dof,x,y,U,xx,yy,Cx,Cy] = solve_lin_rect(N,N,el,el,nh,nh,f,c);

quadplot(x,y,U,xx,yy,Cx,Cy),  show_order(dof,N)

%% Square domain with non-square spectral elements
Nx  = 8;
Ny  = 12;
xel = -1:0.5:1;
yel = [-1,-0.2,0.5,1];
f   = [];
c   = -1;

[dof,x,y,U,xx,yy,Cx,Cy] = solve_lin_rect(Nx,Ny,xel,yel,nh,nh,f,c);

quadplot(x,y,U,xx,yy,Cx,Cy),  show_order(dof,N)

%% Interpolate to a different mesh  (h = 0.1)
h = 0.1;

[Cx,xx] = pbary(Nx,xel,h);  
[Cy,yy] = pbary(Ny,yel,h);  

quadplot(x,y,U,xx,yy,Cx,Cy)

%% Rectangular domain with non-uniform spectral order (Nx,Ny)
Nx  = [ 9,7,10,12 ];
Ny  = 13;
xel = [0, 0.5, 1, 1.8, 3];
yel = [0,1,3,4];
f   = [];
c   = -1;

[dof,x,y,U,xx,yy,Cx,Cy] = solve_lin_rect(Nx,Ny,xel,yel,nh,nh,f,c);

quadplot(x,y,U,xx,yy,Cx,Cy),  show_order(dof,N)

%% Rectangular domain (single element, standard pseudospectral)
Nx  = 25;
Ny  = 20;
xel = [0 3];
yel = [0 4];
f   = [];
c   = -1;

[dof,x,y,U,xx,yy,Cx,Cy] = solve_lin_rect(Nx,Ny,xel,yel,nh,nh,f,c);

quadplot(x,y,U,xx,yy,Cx,Cy),  show_order(dof,N)

%% Square domain, asymmetric forcing (Trefethen, Problem 16)
N  = 8;
el = -1:0.5:1;
nh = 40;
f  = @f16;
c  = 0;

[dof,x,y,U,xx,yy,Cx,Cy] = solve_lin_rect(N,N,el,el,nh,nh,f,c);

quadplot(x,y,U,xx,yy,Cx,Cy),  show_order(dof,N)

%% Laplace equation on square domain with variable boundary-values
% This is Trefethen's Problem 36 (Chapter 13).
N   = 10;
el  = -1:0.5:1;
nh  = 40;
f   = [];
c   = 0;
fbv = @fbv36;

[dof,x,y,U,xx,yy,Cx,Cy] = solve_lin_rect(N,N,el,el,nh,nh,f,c,fbv);

quadplot(x,y,U,xx,yy,Cx,Cy),  show_order(dof,N)

