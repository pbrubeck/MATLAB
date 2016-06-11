function w = chebfftD(u, dim)
% Calculates the partial derivative of u along the dimension dim of 
% a Chebyshev grid. Becomes undefined at the boundary.
N=size(u, dim)-1;
th=(1:N-1)'*pi/N;
D=1i*[0:N-1, 0, 1-N:-1]';
if(dim>1)
    th=reshape(th, [ones(1, dim-1), N-1]);
    D=reshape(D, [ones(1, dim-1), 2*N]);
end
mid=repmat({':'}, 1, ndims(u)); mid{dim}=2:N;
v_hat=fft(cat(dim, u, flip(u(mid{:}), dim)), [], dim);
W=ifft(bsxfun(@times, D, v_hat), [], dim);
w=zeros(size(u));
w(mid{:})=bsxfun(@times, -1./sin(th), W(mid{:}));
end