function uT=unitTangent(v)
% Calculates the unit vectors of the spectral derivate of a curve.
uT=normr(spPartialD(v));
end