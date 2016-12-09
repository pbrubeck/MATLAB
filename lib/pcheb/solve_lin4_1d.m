function  [x,u,xx,uu] = solve_lin4_1d(N,xel,xnh,func)
%solve_lin4_1d  Solve the 1D biharmonic eqn over the interval [-1,1]:
%                   Uxxxx = func(x)  in |x| < 1,
%                       U = 0        at |x| = 1,
%                      Ux = 0        at |x| = 1.

% choose between Chebyshev (0) and Legendre (1) points
xflag = 0;

% choose symmetry of D-matrices - any one of {0,-1,+1}
dflag = 0;

% must use 4th-order mortar BCs
mflag = 1;

[x,xx,C,~,D2,m,M,D3,D4] = pcheb(N,xel,xnh,xflag,dflag,mflag);

n = length(x);  i = 2:n-1;

% Implement Neumman BCs as follows (Trefethen, Chapter 14):
% Let 
%      U(x) = (1 - x^2) p(x).
% By enforcing
%        V=0 at |x|=1,
% we simultaneously enforce Ux=0 at |x|=1.
% Note that
%       Uxxxx = (1 - x^2)*Vxxxx - 8*x*Vxxx - 12*Vxx.

% Mapping from U to V
S = diag([ 0; 1./(1 - x(i).^2); 0; ]); 

DD4 = (diag(1-x.^2)*D4 - 8*diag(x)*D3 - 12*D2)*S;

% 1D Laplace equation
A = DD4;  b = func(x);  

% mortar BCs - continuity of first derivative between elements
A(m,:) = 0;  A = A + M;  b(m) = 0;

% solve and restore boundary points
u = [ 0; A(i,i)\b(i); 0 ];   

% interpolate to fine grid
uu = C*u;

end
