function grad=fftGrad(u)
% Calculates the gradient if u is a scalar field, or the Jacobian
% matrix if u is a vector field. Must be scaled by 2pi/L.
grad=[];
for n=1:ndims(u)
    grad=cat(ndims(u)+1, grad, fftD(u, n));
end
end