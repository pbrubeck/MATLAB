function w = chebfftD2(u, dim)
% Calculates the second partial derivative of u along the dimension dim of 
% a Chebyshev grid.
N=size(u, dim)-1;
th=(1:N-1)'*pi/N;
D=1i*[0:N-1, 0, 1-N:-1]';
if(dim>1)
    D=reshape(D, [ones(1, dim-1), 2*N]);
    th=reshape(th, [ones(1, dim-1), N-1]);
end
index=repmat({':'}, 1, ndims(u)); index{dim}=2:N;
v=cat(dim, u, flip(u(index{:}), dim));
v_hat=fft(v, [], dim);
W1=ifft(bsxfun(@times, D, v_hat), [], dim);
W2=ifft(bsxfun(@times, D.^2, v_hat), [], dim);
w=zeros(size(u));
c=cos(th); s=sin(th);
w(index{:})=bsxfun(@times, s.^-2, W2(index{:}))-bsxfun(@times, c.*s.^-3, W1(index{:}));
end