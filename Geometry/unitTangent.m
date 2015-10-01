function uT=unitTangent(v)
% Calculates the unit vectors of the spectral derivate of a curve.
w=spPartialD(v);
uT=bsxfun(@rdivide, w, sqrt(dot(w, w, ndims(v))));
end