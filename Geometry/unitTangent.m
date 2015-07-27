function uT=unitTangent(v)
% Calculates the unit vectors of the spectral derivate of a curve.
w=spectralD(v,1,1);
uT=bsxfun(@ldivide, sqrt(sum(w.*conj(w), ndims(v))), w);
end