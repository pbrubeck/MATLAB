%% Program 38 of Trefethen (Chapter 14) - ODE with analytic solution
  
% specify spectral order of each element
  N = 10; 
  
% specify domain and element boundaries
  xel = -1:0.5:1;  
  
% force-function
  func = inline('exp(x)');

  xnh = 40; 

  [x,u,xx,uu] = solve_lin4_1d(N,xel,xnh,func);

% plot solution twice, with markers at collocation points
  plot(x,u,'b.',xx,uu,'b-'),  grid on
