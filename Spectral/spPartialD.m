function w=spPartialD(u, dim, n)
% Calculates the spectral n-th partial derivate of a peridic sample along 
% the dimension dim. Must to be scaled by (2pi/L)^n.
M=size(u,dim);
D=1i*[0:M/2-1 -M/2:-1]';
if(dim>1)
    D=reshape(D, [ones(1, dim-1), M]);
end
u_hat=fft(u, [], dim);
w=ifft(bsxfun(@times, u_hat, D.^n), [], dim);
if(isreal(u))
    w=real(w);
end
end