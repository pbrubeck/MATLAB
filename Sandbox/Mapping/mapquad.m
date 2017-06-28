function [f, J, G11, G12, G22] = mapquad(Z, iflag)
% Mapping from [-1,1]^2 to any quadrilateral
if nargin<2
    iflag='noinv';
end

Z=reshape(Z,[2,2]);
T=[1, 1; 1,-1]/2;
A=T'*Z*T;
B=imag(conj(A(:,2))*A(2,:));
E=real(A(2,:)'*A(2,:));
F=real(A(:,2)*A(:,2)');
G=real(conj(A(:,2))*A(2,:));

f=@(x,y) A(1,1)+A(2,1)*x+A(1,2)*y+A(2,2)*(x.*y);
J=@(x,y) B(1,1)+B(2,1)*x+B(1,2)*y+B(2,2)*(x.*y);

if nargout>2
    if strcmp(iflag,'inv')
        G11=@(x,y) (F(1,1)+(F(2,1)+F(1,2))*x+F(2,2)*x.^2)./J(x,y).^2;
        G22=@(x,y) (E(1,1)+(E(2,1)+E(1,2))*y+E(2,2)*y.^2)./J(x,y).^2;
        G12=@(x,y) -(G(1,1)+G(2,1)*x+G(1,2)*y+G(2,2)*(x.*y))./J(x,y).^2;
    else
        G11=@(x,y) E(1,1)+(E(2,1)+E(1,2))*y+E(2,2)*y.^2;
        G22=@(x,y) F(1,1)+(F(2,1)+F(1,2))*x+F(2,2)*x.^2;
        G12=@(x,y) G(1,1)+G(2,1)*x+G(1,2)*y+G(2,2)*(x.*y);
    end
end
end