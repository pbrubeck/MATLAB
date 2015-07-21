function S=tubular(C, r, m)
% Generates a tubular surface from a curve using the Frenet-Serret frame.
uT=unitTangent(C);
uN=unitTangent(uT);
uB=bsxfun(@cross, uT, uN);

v=2*pi*(0:m-1)'/m;
S=kron(C, ones(m,1))+kron(uN, r*cos(v))+kron(uB, r*sin(v));
S=reshape(S, [m, size(C)]);
end