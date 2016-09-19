 function dw = sgrhs(t,w,D)

%  The function dw = sgrhs(t,w,D) computes the right-hand side
%  of the sine-Gordon system with the aid of an NxN differentiation matrix D.
 
%  J.A.C. Weideman, S.C. Reddy 1998
%  Modified for MATLAB R13 by JACW in April 2003. 

 N = length(w)/2;
 u = w(1:N); v = w(N+1:2*N);

du = v;
dv = D*u-sin(u);

dw = [du; dv];




