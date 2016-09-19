
%  The script file matplot.m plots the characteristic curves of
%  Mathieu's equation using a 32x32 Fourier differentiation matrix.

%  J.A.C. Weideman, S.C. Reddy 1998

[x, D] = fourdif(32,2); 
for q  = 0.1:0.1:12                         % For loop over q values
     a = eig(2*q*diag(cos(2*x))-D);         % Compute eigenvalues (period 2 pi)
plot(q+i*a,'o','MarkerSize',2); hold on;
end
axis([0 12 -10 30])                         % Zoom in

xlabel('q','FontSize',16)
ylabel('a','FontSize',16)
