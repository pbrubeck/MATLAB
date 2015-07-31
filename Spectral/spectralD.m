function w = spectralD(u, n, dim)
% Calculates the spectral n-th derivate of a peridic sample along the
% dimension dim. Actual derivate has to be scaled by 2pi/L.
if(nargin==1)
    n=1; dim=1;
end
if(nargin==2)
    dim=1;
end
N=size(u, dim);
k=1i*[0:N/2-1, 0, -N/2+1:-1]';
if(dim>1)
    k=reshape(k, [ones(dim-1), N]);
end
w=ifft(bsxfun(@times, k.^n, fft(u, [], dim)), [], dim);
end