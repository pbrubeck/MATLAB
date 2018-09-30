function [S1,S2,S3,S4] = stokesparams(Ex,Ey)
S1=abs(Ex).^2+abs(Ey).^2;
S2=abs(Ex).^2-abs(Ey).^2;
S3=2*real(conj(Ex).*Ey);
S4=2*imag(conj(Ex).*Ey);
end