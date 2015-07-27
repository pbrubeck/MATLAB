function uN=surfNormal(S)
% Calcules the normal vector field of a given surface through spectral
% differentiation

Su=spectralD(S,1,1);
Sv=spectralD(S,1,2);
N=cross(Su, Sv);
dA=sqrt(sum(N.*conj(N), ndims(N)));

uN=bsxfun(@ldivide, dA, N);
end