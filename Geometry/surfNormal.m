function uN=surfNormal(S)
% Calcules the normal vector field of a given surface through spectral
% differentiation
n=size(S,1);
S_hat=fft(S);
Su_hat=bsxfun(@times, 1i*[0:n/2-1, 0, -n/2+1:-1]', S_hat);
Su=real(ifft(Su_hat));

m=size(S,2);
S_hat=fft(permute(S, [2 1 3]));
Sv_hat=bsxfun(@times, 1i*[0:m/2-1, 0, -m/2+1:-1]', S_hat);
Sv=permute(real(ifft(Sv_hat)), [2 1 3]);

uN=bsxfun(@cross, Su, Sv);
uN=bsxfun(@ldivide, sqrt(sum(uN.*conj(uN), ndims(uN))), uN);
end

