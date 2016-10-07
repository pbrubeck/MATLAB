function [] = sqhole( N )
% Conformal mapping square with circular hole
a=0.5; b=1;
[Dr, r]=chebD(N); r=a+(b-a)*(r-1)/2; Dr=2/(b-a)*Dr; Drr=Dr*Dr;
[Dtt, t]=fourD2(N);

f=mod(t/(pi/4)+1,2)-1;
wb=[b*1i.^round(4*t/(2*pi)).*(1+1i*2/pi*asin(f)); a*exp(1i*t)];
RHS=-Drr(:,[1 end])*wb;

% Solve Laplace for the conformal mapping
w=zeros(size(Drr)); w([1 end],:)=wb;
w(2:end-1,:)=sylvester(Drr(2:end-1,2:end-1), Dtt', RHS(2:end-1,:));

% Compute jacobian determinant
J=abs(Dr*w).^2;

% Solve Poisson 
bc=zeros(2, N);
F=-ones(N);
RHS=J.*F-Drr(:,[1 end])*bc;
psi=zeros(size(Drr)); 
psi([1 end],:)=bc;
psi(2:end-1,:)=sylvester(Drr(2:end-1,2:end-1), Dtt', RHS(2:end-1,:));

u=real(w);
v=imag(w);
surf(u(:,[1:end,1]), v(:,[1:end,1]), psi(:,[1:end,1]));
axis square; shading interp; colormap(jet(256));
end

