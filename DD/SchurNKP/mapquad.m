function [zz,jac,g11,g12,g22]=mapquad(Z,x,y)
% Bilinear map from [-1,1]^2 to arbitrary quadrilaterals
% Z=[z(NE), z(NW); z(SE), z(SW)]
T=[1,1;1,-1]/2;
A=T'*reshape(Z,[2,2])*T;
J=imag(conj(A(:,2))*A(2,:));
E=real(A(:,2)*A(:,2)');
F=real(conj(A(:,2))*A(2,:));
G=real(A(2,:)'*A(2,:));

% Transformation
zz=A(1,1)+A(2,1)*x+A(1,2)*y+A(2,2)*x.*y;

% Jacobian
jac=J(1,1)+J(2,1)*x+J(1,2)*y;

% Metric
g11=E(1,1)+(E(2,1)+E(1,2))*y+E(2,2)*y.^2;
g12=F(1,1)+F(1,2)*x+F(2,1)*y+F(2,2)*x.*y;
g22=G(1,1)+(G(2,1)+G(1,2))*x+G(2,2)*x.^2;
end