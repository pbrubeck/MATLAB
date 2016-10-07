function [ ] = sq2circ(n)
% Conformal mapping from [-1,1]^2 to the unit circle
[D,x]=chebD(n); D2=D*D;
[xx,yy]=ndgrid(x);

% Boundary conditions
th=-pi/4*x';
bx=[exp(1i*th); exp(1i*(pi-th))];
by=[exp(1i*(3*pi/2-th')), exp(1i*(pi/2+th'))];
RHS=-D2(:,[1 end])*bx-by*D2(:,[1 end])';

% Solve Laplace
w=zeros(size(D2)); w([1 end],:)=bx; w(:,[1 end])=by;
w(2:end-1,2:end-1)=sylvester(D2(2:end-1,2:end-1), D2(2:end-1,2:end-1)', RHS(2:end-1,2:end-1));

im=imread('airplane.tiff');
figure(1);
imshow(im);

[xq,yq]=ndgrid(linspace(-1,1,size(im,1)),linspace(-1,1,size(im,2)));
wq=interp2(xx',yy',w.',xq,yq,'spline').';
figure(2);
h=surf(real(wq),imag(wq),zeros(size(wq)));
set(h,'CData',im); shading interp; colormap(gray);
axis square; view(0,90);
end

