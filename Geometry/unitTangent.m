function uT=unitTangent(v)
% Calculates the unit vectors of the spectral derivate of a curve.
uT=normr(fftD(v,1,1));
end