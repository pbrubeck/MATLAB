function div = spDivergence(u)
% Calculates the spectral divergence of a vector field, actual divergence
% must be scaled by 2pi/L.
div=spPartialD(u(:,:,:,1),1)+spPartialD(u(:,:,:,2),2)+spPartialD(u(:,:,:,3),3);
end