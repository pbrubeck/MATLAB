function div = spDivergence(u)
% Calculates the spectral divergence of a vector field, actual divergence
% must be scaled by 2pi/L.
N=size(u);
i=[0:N(1)/2-1, 0, -N(1)/2+1:-1];
j=[0:N(2)/2-1, 0, -N(2)/2+1:-1];
k=[0:N(3)/2-1, 0, -N(3)/2+1:-1];
[ii,jj,kk]=meshgrid(i, j, k);
omega=cat(4, ii, jj, kk);
div=spDot(1i*omega, u);
end