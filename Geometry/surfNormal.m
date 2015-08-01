function uN=surfNormal(S)
% Calcules the normal vector field of a given surface through spectral
% differentiation
N=cross(spPartialD(S,1), spPartialD(S,2));
dA=sqrt(dot(N, N, 3));
uN=bsxfun(@rdivide, N, dA);
end