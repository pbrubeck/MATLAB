function w = chebfftD2(u, dim)
% Calculates the second partial derivative of u along the dimension dim of 
% a Chebyshev grid. Becomes undefined at the boundary.
N=size(u, dim)-1;
th=(1:N-1)'*pi/N;
D=1i*[0:N, 1-N:-1]';
if(dim>1)
    D=reshape(D, [ones(1, dim-1), 2*N]);
    th=reshape(th, [ones(1, dim-1), N-1]);
end
mid=repmat({':'}, 1, ndims(u)); mid{dim}=2:N;
v_hat=fft(cat(dim, u, flip(u(mid{:}), dim)), [], dim);
W1=ifft(bsxfun(@times, D, v_hat), [], dim);
W2=ifft(bsxfun(@times, D.^2, v_hat), [], dim);
w=zeros(size(u));
c=cos(th); s=sin(th);
w(mid{:})=bsxfun(@times, s.^-2, W2(mid{:}))-bsxfun(@times, c.*s.^-3, W1(mid{:}));
end