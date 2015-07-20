function uT=unitTangent(v)
% Calculates the unit vectors of the spectral derivate of a curve.
n=size(v,1);
v_hat=fft(v);
w_hat=bsxfun(@times, 1i*[0:n/2-1, 0, -n/2+1:-1]', v_hat);
w=real(ifft(w_hat));
uT=w./repmat(sqrt(sum(w.*conj(w),2)),[1 3]);
end