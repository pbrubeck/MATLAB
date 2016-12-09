function [D,x] = leg(N,xmax,xmin)
%LEG  Legendre spectral collocation scheme in one dimension.
%       [Adapted from "legDc.m" contributed to Matlab Central
%        by Greg von Winckel, gregvw@chtm.unm.edu]
% 
% [D,x] = leg(N,xmax,xmin)
%
%  x - Legendre collocation points in ascending order over [xmin,xmax]
%      [column-vector of length N+1]
%
% D - differentiation matrix
%      [square matrix, size N+1]
%
% See also LEGDC, MYCHEB, PCHEB, DIFFBARY.


    M=N+1;

    [~,xc] = cheb(N);

    % Uniform nodes
    xu = linspace(-1,1,M)';

    % Make a close first guess to reduce iterations
    if N < 3
        x = xc;
    else
        x = xc + sin(pi*xu)./(4*N);
    end

    % The Legendre Vandermonde Matrix
    P = zeros(M,M);

    % Compute P_(N) using the recursion relation
    % Compute its first and second derivatives and 
    % update x using the Newton-Raphson method.

    xold=2;
    
    % collocation points on [-1,1] in reverse order +1 ==> -1
    while max(abs(x-xold))>eps

        xold=x;

        P(:,1)=1;    P(:,2)=x;

        for k=2:N
            P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
        end

        x = xold-( x.*P(:,M)-P(:,N) )./( M*P(:,M) );
        
    end
    
    % collocation points on [xmin,xmax] in forward order xmin ==> xmax
    x = (xmin+xmax)/2 - (xmax - xmin)*x/2;
    
    D = diffbary(x);


end

