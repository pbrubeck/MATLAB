function [m,M,xv,yv] = mort_rect(x,y,mx,my,Mx,My)
%MORT_RECT  Mortar BCs for spectral elements over a rectangular domain.
%
% [m,M,xv,yv] = mort_rect(x,y,mx,my,Mx,My)
%
%  x - x-wise collocation points (length nx):
%
%         [x,~,~,~,~,mx,Mx] = pcheb(Nx,xel);
%
%  y - y-wise collocation points (length ny):
%
%         [y,~,~,~,~,my,My] = pcheb(Ny,yel);
%
% xv, yv - collocation points as column vectors of length nxy = nx*ny
%
%      m - mortar-points, as sparse logical column vector of length nxy
%
%      M - mortar BCs, as sparse square matrix of dimension nxy
%
% See also PCHEB, PBARY, MORT_SQ.



    nx = length(x);  Ix = speye(nx);
    ny = length(y);  Iy = speye(ny);

    [X,Y] = meshgrid(x,y);  xv = X(:);  yv = Y(:);
    
    b = (xv==x(1) | xv==x(end) | yv==y(1) | yv==y(end)); % boundary
    
    M = kron(Mx,Iy) + kron(Ix,My);

    m  = (ismember(xv,x(mx)) | ismember(yv,y(my)));
    
    m = m & ~b;   M(b,:) = 0;  % remove mortar BCs on boundary

end
