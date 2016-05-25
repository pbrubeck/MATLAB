function div = fftDiv(u)
% Calculates the spectral divergence of a vector field, actual divergence
% must be scaled by 2pi/L.
div=fftD(u(:,:,:,1),1)+fftD(u(:,:,:,2),2)+fftD(u(:,:,:,3),3);
end