function grad=spGradient(Phi)
% Calculates the spectral divergence of a vector field, actual divergence
% must be scaled by 2pi/L.
N=size(Phi);
i=[0:N(1)/2-1, 0, -N(1)/2+1:-1];
j=[0:N(2)/2-1, 0, -N(2)/2+1:-1];
k=[0:N(3)/2-1, 0, -N(3)/2+1:-1];
[ii,jj,kk]=meshgrid(i, j, k);
Phi_hat=fftn(Phi);
grad(:,:,:,1)=1i*ifftn(ii.*Phi_hat);
grad(:,:,:,2)=1i*ifftn(jj.*Phi_hat);
grad(:,:,:,3)=1i*ifftn(kk.*Phi_hat);
end