function [ ] = pwePML2( N )
% Solves the paraxial wave equation in 2D using a Perfectly Matched Layer for BCs.
% Uses Hermite spectral methods in two spatial dimensions
% and classic Runge-Kutta for z evolution of second order PDE.
lambda=100;
k=2*pi/lambda;
[D, x]=hermD(N);

[xx,yy]=ndgrid(x);
dz=6/N^2;

% Layer
xl=3/4*x(end);
layer=abs(x)>xl;
roi=~layer;

sig=zeros(N,1);
sig(layer)=2/dz*((abs(x(layer))-xl)/(x(end)-xl)).^3;

% Initial conditions
u0=(cos(2*pi*(xx+yy))+cos(2*pi*(xx-yy))).*exp(-(xx.^2+yy.^2)/2);
vx=zeros(N);
vy=zeros(N);
psi=zeros(N);
w=cat(3,u0,vx,vy,psi);

function wt=partialTime(w)
wt=w;
u=w(:,:,1);
vx=D*u -dgmm(sig , w(:,:,2));
vy=u*D'-dgmm(sig', w(:,:,3));
wt(:,:,2)=vx;
wt(:,:,3)=vy;
wt(:,:,4)=-dgmm(sig,dgmm(sig',u))+dgmm(sig',vx)+dgmm(sig,vy);
wt(:,:,1)=1i/(2*k)*(D*vx+vy*D')-dgmm(sig,u)-dgmm(sig',u)+w(:,:,4);
end

function w=solveRK4(w, dz)
k1=dz*partialTime(w);
k2=dz*partialTime(w+k1/2);
k3=dz*partialTime(w+k2/2);
k4=dz*partialTime(w+k3);
w=w+(k1+2*k2+2*k3+k4)/6;
end

figure(1);
h=image(x(roi), x(roi), abs(w(roi,roi,1)).^2);
map=hsv(256); colormap(map); shading interp;
view(2); axis square; 

nframes=10000;
for i=1:nframes
    w=solveRK4(w,dz);
    w([1 end],:,:)=0;
    w(:,[1 end],:)=0;
    if(mod(i,5)==1)
        E=w(roi,roi,1);
        
        V=mat2gray(abs(E).^2);
        V=cat(3,V,V,V);
        phase=angle(E)/(2*pi);
        phase=phase-floor(phase);
        H=ind2rgb(uint8((length(map)-1)*phase), map);
        
        h.CData=H.*V;
        drawnow;
    end
end
end

function A=dgmm(x,A)
A=bsxfun(@times, x, A);
end