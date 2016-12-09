function D = diffbary(x)
%DIFFBARY  Spectral differentation matrix for given collocation points.
%
% See also BARYCENTRIC.


    x   = x(:);
    M   = length(x);   % M = N+1, where N is the polynomial order
    c   = zeros(M,1);
    X   = repmat(x,1,M);
    dX  = X - X';
    Eps = eye(M);  % avoid division-by-zero
    
    for n = 1:M
        xx   = x([ 1:n-1, n+1:M ]);  
        c(n) = prod(x(n) - xx);
    end
    
    % set off-diagonal entries
    D = (c * (1./c)')./(dX + Eps);
    
    % set on-diagonal entries
    %
    % METHOD 1 [Direct, formula-based]:
    %    D(n,n) = sum(1./xx);  [using XX as defined above]
    %
    % METHOD 2 [Indirect, ensures all row-sums are zero]:
         D = D - diag(sum(D,2));

end

