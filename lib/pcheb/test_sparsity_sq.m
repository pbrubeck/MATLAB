%% Single square element
N   = 20;
xel = [0,1]; % any ordered vector of length M2 will give the same result
xnh = 50;    % interpolation vector, determines number of C-rows

[x,~,~,D,D2] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

Dxx = kron(D2, I);
Dyy = kron( I,D2);
Dxy = kron( D, D);
%
% Equivalently:
%   Dxy = kron(D,I)*kron(I,D)
%   Dxy = kron(I,D)*kron(D,I)
%
% In general,
%   kron(A*C,B*D) = kron(A,B)*kron(C,D)
% for any quartet [A,B,C,D] of square or (compatible) rectangular matrices.
%
% This is a generic property of the matrix Kronecker product:
%
%   http://en.wikipedia.org/wiki/Kronecker_product

str = '  [M = 1]';

A=Dxx; 
subplot(2,3,1),  spy(A),  d = spdensity(A);
title(['Dx, Dxx',str,'  --  density = ',num2str(d*N,'%5.3f'),' / N'])

A=Dxy;
subplot(2,3,2),  spy(A),  
title(['Dxy',str,'  --  density = 1'])

A=Dyy;
subplot(2,3,3),  spy(A),  
d = spdensity(A);
title(['Dy, Dyy',str,'  --  density = ',num2str(d*N,'%5.3f'),' / N'])

%% Mesh of M x M elements
N   = 5;
M   = 4;    % number of elements
xel = 0:M;  % any ordered vector of length M+1 will give the same result
xnh = 50;   % interpolation vector, determines number of C-rows

str = ['  [M = ',int2str(M),']'];

Msq = M^2;

[x,~,~,D,D2] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

Dxx = kron(D2, I);
Dyy = kron( I,D2);
Dxy = kron( D, D);


A=Dxx; 
subplot(2,3,4),  spy(A),  
d = spdensity(A);
title(['Dx, Dxx',str,'  --  density = ',num2str(d*Msq*N,'%5.3f'),' / M^2N'])

A=Dxy;
subplot(2,3,5),  spy(A),  
d = spdensity(A);
title(['Dxy',str,'  --  density = ',num2str(d*Msq,'%5.3f'),' / M^2'])

A=Dyy;
subplot(2,3,6),  spy(A),  
d = spdensity(A);
title(['Dy, Dyy',str,'  --  density = ',num2str(d*Msq*N,'%5.3f'),' / M^2 N'])

