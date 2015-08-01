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
Dn=(1i*[0:N/2-1, 0, -N/2+1:-1]').^n;
if(dim>1)
    Dn=reshape(Dn, [ones(1, dim-1), N]);
end
w=ifft(bsxfun(@times, Dn, fft(u, [], dim)), [], dim);
end