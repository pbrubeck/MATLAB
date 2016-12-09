function [x,xx,C,D1,D2,mort,Mort,D3,D4] = ...
                                 pcheb(Nvec,xel,xnh,xflag,dflag,mflag)
% PCHEB  Spectral-collocation hp-elements in one dimension:
%
%      p-element method ==> vary Nvec to improve resolution
%      h-element method ==> vary  XEL to improve resolution
%     hp-element method ==> vary Nvec and XEL to improve resolution
%
% [x,xx,C,D1,D2,mort,Mort,D3,D4] = pcheb(Nvec,xel,xnh,xflag,dflag,mflag)
%
%
% == INPUTS [mandatory] ==
%
%    xel - vector of length M+1 specifying boundaries of M elements
%
%   Nvec - polynomial order, in one of the following formats:
%           - integer N>0  (applied to all elements)
%           - integer-valued vector Nvec>0 of length M.
%
% == OUTPUT ==
%
%    x - collocation points
%        [column-vector of length 'n', n = M*N + 1]
%
%   xx - interpolation points 
%        [column-vector of length 'm']
%
%    C - interpolation matrix X => XX  
%        [size m-by-n, sparse]
%
% [D1,D2,D3,D4] - differentation matrices of order 1,2,3,4 respectively
%        [size n-by-n, sparse]
%
% Mort - mortar boundary conditions (BCs)
%        [size n-by-n, sparse]
%
% mort - mortar-points, i.e. collocation points at element boundaries
%        [logical column-vector, length 'n']
%
%
% == INPUTS [optional] ==
%
%   xnh - XX, specified in one of the following formats:
%              [empty vector] - default
%        [vector, length 'm'] - prescribed
%         [scalar, 0 < h < 1] - grid spacing (equispaced)
%           [integer, n >= 1] - n+1 points, equispaced
%
% xflag - [binary] choice of collocation points:
%                     0/false ==>  Chebyshev
%                     1/true  ==>  Legendre
%
% dflag - [ternary] treatment of mortar-points in [D1,D2,D3,D4]:
%           +1 ==>  evaluate in right-hand element 
%           -1 ==>  evaluate in  left-hand element 
%            0 ==>  average of left- and right-hand derivatives
%
% mflag - [binary] regularity of mortar BCs
%                       0/false ==>  2nd-order [1 continuous derivative ]
%                       1/true  ==>  4th-order [3 continuous derivatives]
%
% NOTE:
%   DFLAG is usually of no practical significance, since mortar-values of 
%     differentiation matrices are typically over-written by mortar BCs.
%   The 3 cases [-1,0,+1] are then equivalent.
%   Occasionally, however, the correct choice of DFLAG is crucial - see,
%     for example, SOLVE_LIN_LWISE_BC2 and SOLVE_LIN_SQ_LREFINE.
%
% See also MYCHEB, LEG, BARYCENTRIC, PBARY, MORT_RECT, MORT_SQ,
%          SOLVE_LIN_LWISE_BC2 and SOLVE_LIN_SQ_LREFINE.


if nargin < 3 || isempty(xnh),  xnh = 50;  end 

if nargin < 4,  xflag = 0;  end
if nargin < 5,  dflag = 0;  end
if nargin < 6,  mflag = 0;  end


M = length(xel) - 1;

if M < 1,
    error('Arg #2 (''xel'') should be a vector of length >= 2.')
end

Mn = length(Nvec);

if Mn==1
    Nvec = Nvec*ones(1,M);
elseif Mn < M
    error(['Arg #1 should be a scalar or a ',int2str(M),'-vector.'])
end

xel = sort(xel);  x0 = xel(1);  x1 = xel(end);

if isscalar(xnh)

    if xnh <= 0,  error('Arg 3 faulty.'),  end

    if xnh >= 1,
        n  = round(xnh);
        xx = x0 + (x1 - x0)*(0:n)'/n;
    else
        h  = xnh;
        xx = (x0:h:x1)';
    end

else
    xx = sort(xnh(:));
end
    
%   i - index of X with    duplicate mortar-points  [length n0]
%   j - index of X without duplicate mortar-points  [length n ]
% i2j - mapping from i to j                         [length n0]
    
    Zero = sparse(0,0);
        
    x = [];  j=1;  i2j = [];
    
    C = Zero;  D1 = Zero;  D2 = Zero;  D3 = Zero;  D4 = Zero;
    
    for m = 1:M,
        
        xm0 = xel(m);  xm1 = xel(m+1);  N = Nvec(m);
        
        i2j = [ i2j  j+(0:N) ];   j=j+N;
        
        if xflag,
            [Dm,xm] =    leg(N,xm1,xm0); 
        else
            [Dm,xm] = mycheb(N,xm1,xm0); 
        end
        
        if m==1
            xxm = xx(xx >= xm0  &  xx <= xm1);
        else
            xxm = xx(xx >  xm0  &  xx <= xm1);
        end
        
        Cm = barycentric(xm,xxm);
        
      % accumulate results
        x  = [ x; xm ];
        C  = blkdiag( C,Cm  );
        D1 = blkdiag(D1,Dm  );
        D2 = blkdiag(D2,Dm^2);
        D3 = blkdiag(D3,Dm^3);
        D4 = blkdiag(D4,Dm^4);
         
    end
    
    n0   = sum(Nvec+1);
    n    = sum(Nvec)+1;
    ii   = 1:n0;
    one  =  ones(1,n0);
    zero = zeros(1,n0);
    NN   = Nvec(1:end-1);
    jm   = cumsum(NN)+1; %     unique mortar-indices
    iml  = cumsum(NN+1); %  left-hand mortar-indices
    imr  = iml+1;        % right-hand mortar-indices
    
    % remove duplicate mortar-columns [additive merge]
    kk = one;
    S  = sparse(i2j,ii,kk);  % N-by-N0 matrix [N0 non-zero entries]
    C  =  C*S';
    D1 = D1*S';
    D2 = D2*S';
    D3 = D3*S';
    D4 = D4*S';
    
    % mortar boundary conditions [subtractive row-wise merge of D1]
    kk = zero;  kk(iml) = +1;  kk(imr) = -1;  
         % OR:  kk(iml) = -1;  kk(imr) = +1;
         
    S    = sparse(i2j,ii,kk);  % N-by-N0 matrix
    Mort = S*D1;               % ensure D1-continuity at mortars
    mort = false(n,1);   mort(jm) = true;
    
    if mflag,
        
        Mort2 = S*D2;  % ensure D2-continuity at mortars
        Mort3 = S*D3;  % ensure D3-continuity at mortars
        
        jm2 = jm-1;  % place D2-continuity eqn adjacent to mortars (LHS)
        jm3 = jm+1;  % place D3-continuity eqn adjacent to mortars (RHS)
        
        Mort(jm2,:) = Mort2(jm,:);   mort(jm2) = true;
        Mort(jm3,:) = Mort3(jm,:);   mort(jm3) = true;
        
    end
    
    % remove duplicate mortar-rows [semi-additive merge]
    s  = sign(dflag);
    kk = one;   kk(iml) = (1-s)/2;  kk(imr) = (1+s)/2;
    S  = sparse(i2j,ii,kk);  % N-by-N0 matrix
    D1 = S*D1; 
    D2 = S*D2;
    D3 = S*D3; 
    D4 = S*D4;
    
    x = S*x;   % remove duplicate mortar-points

end

