function w=spPartialD(u, dim, n)
% Calculates the spectral n-th partial derivate of a peridic sample along 
% the dimension dim. Must to be scaled by 2pi/L.
if(nargin==1)
    n=1; dim=1;
end
if(nargin==2)
    n=1;
end
N=size(u, dim);
D=1i*[0:N/2-1, 0, -N/2+1:-1]';
if(dim>1)
    D=reshape(D, [ones(1, dim-1), N]);
end
w=ifft(bsxfun(@times, D.^n, fft(u, [], dim)), [], dim);
end