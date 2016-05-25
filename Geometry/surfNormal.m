function uN=surfNormal(S)
% Calcules the normal vector field of a given surface through spectral
% differentiation
N=cross(fftD(S,1,1), fftD(S,2,1));
dA=sqrt(dot(N, N, 3));
uN=bsxfun(@rdivide, N, dA);
end