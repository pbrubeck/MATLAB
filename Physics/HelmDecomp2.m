function [Phi, A] = HelmDecomp2(N)
% Helmholtz decomposition of a 2D field

N([1 2])=N;
[Dx,x]=chebD(N(1)); Dxx=Dx*Dx;
[Dy,y]=chebD(N(2)); Dyy=Dy*Dy;
[xx,yy]=ndgrid(x,y);

zz=xx+1i*yy;
F=-1i*(zz.^4).*exp(-10*abs(zz).^2);
P=zeros(N);
P(2:end-1,2:end-1)=sylvester(Dxx(2:end-1,2:end-1), Dyy(2:end-1,2:end-1)', -F(2:end-1,2:end-1));

Phi=Dx*real(P)+imag(P)*Dy'; % Phi=div(P)
A=Dx*imag(P)-real(P)*Dy';   % A=curl(P)=div(-1i*P)

% F=-grad(Phi)+curl(A)=grad(-Phi-1i*A)
% grad(f)=(Dx*f)+1i*(f*Dy')
end