 %% Program 13 of Trefethen (Chapter 7) - ODE with analytic solution
 
% specify spectral order of each element
  N = 8; 
  
% specify domain and element boundaries
  xel = -1:0.5:1; 
  
% force-function
  func = inline('exp(4*x)');  % Problem 13, Trefethen Chapter 7

% specify output grid - in any one of 4 ways
  xnh = 0.05;       % grid-spacing

 [x,u,xx,uu] = solve_lin_1d(N,xel,xnh,func);

% known exact solution for Problem 13
  uex = (exp(4*xx) - xx.*sinh(4) - cosh(4))/16; 

% plot solution twice, with markers at collocation points
  plot(x,u,'b.',xx,uu,'b-'),  grid on

% check accuracy of solution along output grid
  err = norm(uu - uex, inf);

  title(['max-error = ',num2str(err)])
 
%% Repeat with variable elements
% specify spectral order of each element
  N = [10,8,9]; 
  
% specify domain and element boundaries
  xel = [ -1, -0.2, 0.4, 1 ];  % using 3 elements on [-1,1]
  
% force-function
  func = inline('exp(4*x)');  % Problem 13, Trefethen Chapter 7

 xnh = 50;       % points

 [x,u,xx,uu] = solve_lin_1d(N,xel,xnh,func);

% known exact solution for Problem 13
  uex = (exp(4*xx) - xx.*sinh(4) - cosh(4))/16; 

% plot solution twice, with markers at collocation points
  plot(x,u,'b.',xx,uu,'b-'),  grid on

% check accuracy of solution along output grid
  err = norm(uu - uex, inf);

  title(['max-error = ',num2str(err)])
