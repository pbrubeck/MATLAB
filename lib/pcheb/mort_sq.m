function [m,M,xv,yv] = mort_sq(x,mx,Mx)
%MORT_SQ  Mortar BCs for spectral elements over a square domain.
%
% [m,M,xv,yv] = mort_sq(x,mx,Mx)
%
%  x - 1D collocation points (length n):
%
%         [x,~,~,~,~,mx,Mx] = pcheb(Nx,xel);
%
% xv, yv - collocation points as column vectors of length n^2
%
%      m - mortar-points, as sparse logical column vector of length n^2
%
%      M - mortar BCs, as sparse square matrix of dimension n^2
%
% See also MORT_RECT, PCHEB, PBARY.

    [m,M,xv,yv] = mort_rect(x,x,mx,mx,Mx,Mx);

end
