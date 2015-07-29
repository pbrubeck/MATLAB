function div = spDivergence(u)
% Calculates the spectral divergence of a vector field, actual divergence
% must be scaled by 2pi/L.
N=size(u);
i=[0:N(1)/2-1, 0, -N(1)/2+1:-1];
j=[0:N(2)/2-1, 0, -N(2)/2+1:-1];
k=[0:N(3)/2-1, 0, -N(3)/2+1:-1];
[ii,jj,kk]=meshgrid(i, j, k);
div=1i*ifftn(ii.*fftn(u(:,:,:,1))+jj.*fftn(u(:,:,:,2))+kk.*fftn(u(:,:,:,3)));
end