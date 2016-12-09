 %%
  N    = 10; 
  xel  = [ -1, 0, 1 ]; 
  func = inline('exp(u)'); % Problem 14, Trefethen Chapter 7
  xnh  = 0.05;
  tol  = 1e-12;
  
% specifying Problem 13 of Trefethen on domain [-1,1] ...
 [x,u,xc,uc,n] = solve_nlin_1d(N,xel,xnh,func,tol);

% plot solution twice, with markers at collocation points
  plot(xc,uc,'b.',x,u,'b-'),  grid on
  
  title(['Converged after ',int2str(n),' steps'])
