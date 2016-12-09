%% Matrices for 2nd-order ODEs
% The sparsity pattern depends only on the number of elements, not on 
%   their width, uniformity, or type (Chebyshev vs Legendre).
% D1 and D2 have the same sparsity pattern
N   = 5;
M   = 4;    % number of elements
xel = 0:M;  % any ordered vector of length M+1 will give the same result
xnh = 50;   % interpolation vector, determines number of C-rows

% choose between Chebyshev (0) and Legendre (1) points
xflag = 1;

% how to compute derivatives at internal element boundaries (mortars).
% This is usually irrelevant - mortar BCs take precedence.
% DFLAG=0 mimics the effect of mortar BCs.
dflag = 0; 

[~,~,C,D,D2] = pcheb(N,xel,xnh,xflag,dflag);

A=D2;  % try changing this to D
subplot(2,2,1),  spy(A),  
d = spdensity(A);
title(['Order-2 ODE  --  density = ',num2str(d*M,'%5.3f'),' / M'])

A=C;
subplot(2,2,3),  spy(A),  axis square
d = spdensity(A);
title(['Order-2 interp''n  --  density = ',num2str(d*M,'%5.3f'),' / M'])

%% Matrices for 4th-order ODEs
% The sparsity pattern depends only on the number of elements, not their
%   extent or uniformity.
N   = 5;
M   = 4;    % number of elements
xel = 0:M;  % any ordered vector of length M+1 will give the same result
xnh = 50;   % interpolation vector, determines number of C-rows

% choose between Chebyshev (0) and Legendre (1) points
xflag = 0;

% request 4th-order continuity of elements
mflag = 1;

[~,~,C,D,~,~,Cont] = pcheb(N,xel,xnh,xflag,dflag,mflag);

A=D+Cont;  % mimics the effect of mortar BCs
subplot(2,2,2),  spy(A),  
d = spdensity(A);
title(['Order-4 ODE  --  density = ',num2str(d*M,'%5.3f'),' / M'])

A=C;
subplot(2,2,4),  spy(A),  axis square
d = spdensity(A);
title(['Order-4 interp''n  --  density = ',num2str(d*M,'%5.3f'),' / M'])

