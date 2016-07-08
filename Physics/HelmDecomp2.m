function [Phi, A] = HelmDecomp2(N)
% Helmholtz decomposition of a 2D field

N([1 2])=N;
[Dx,x]=chebD(N(1)); Dxx=Dx*Dx;
[Dy,y]=chebD(N(2)); Dyy=Dy*Dy;
[xx,yy]=ndgrid(x,y);

zz=xx+1i*yy;
F=(zz.^8-0.2^4);
P=zeros(N);
P(2:end-1,2:end-1)=sylvester(Dxx(2:end-1,2:end-1), Dyy(2:end-1,2:end-1)', -F(2:end-1,2:end-1));

Phi=Dx*real(P)+imag(P)*Dy'; % Phi=div(P)
A=Dx*imag(P)-real(P)*Dy';   % A=curl(P)=div(-1i*P)

% F=-grad(Phi)+curl(A)=grad(-Phi-1i*A)
% grad(f)=(Dx*f)+1i*(f*Dy')
f=-Phi-1i*A;
E=(Dx*f)+1i*(f*Dy');

figure(1); clf; hold on;
im1=complex2im(E);
surf(xx,yy,-abs(E),im1); shading interp; axis square; view(2);
colormap(hsv(256));
caxis([-pi,pi]);
colorbar('YTick', linspace(-pi, pi, 5), ...
    'YTickLabel', {'-\pi','-\pi/2','0','\pi/2','\pi'});

[xq,yq]=meshgrid(linspace(-0.9,0.9,32));
uu=interp2(xx',yy',E.',xq,yq); uu=uu./abs(uu);
quiver(xq, yq, real(uu), imag(uu), 'w');
hold off;

figure(2);
surfl(xx,yy,Phi,'light'); title('Scalar potential');
colormap(jet(64)); shading interp; axis square; view(2);

figure(3);
surfl(xx,yy,A,'light'); title('Vector potential (z-component)');
colormap(jet(64)); shading interp; axis square; view(2);

end

function im=complex2im(E)
C=log(E);
A=ind2rgb(uint8(128-128*imag(C)/pi), hsv(256));
B=mat2gray(real(C));
B=cat(3, B, B, B);
im=A.*B;
end