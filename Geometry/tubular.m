function [] = tubular(C, r, m)
% Generates a tubular surface from a curve using the Frenet-Serret frame.
uT=unitTangent(C);
uN=unitTangent(uT);
uB=bsxfun(@cross, uT, uN);

v=2*pi*(0:m-1)'/(m-1);
points=repmat(C, [m,1])+kron(r*cos(v), uN)+kron(r*sin(v), uB);
points=reshape(points, [size(C,1),m,3]);
mesh(points(:,:,1), points(:,:,2), points(:,:,3));
axis equal;
end