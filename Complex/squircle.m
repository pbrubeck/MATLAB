function [] = squircle( N )
% Conformal mapping square with circular hole
a=1/2; b=sqrt(2);
[D, x]=chebD(N); r=a+(b-a)*(x+1)/2; D=2/(b-a)*D; A1=(diag(r)*D)^2;
[A2, t]=fourD2(N);

z0=b*exp(2i*pi*(0:7)'/8);
w0=b*z0./(abs(real(z0))+abs(imag(z0)));

g1=(t<pi/2);
g2=(t>=pi/2 & t<pi);
g3=(t>=pi & t<3*pi/2);
g4=(t>=3*pi/2);

H1=mobius(z0(1:3), w0(1:3));
H2=mobius(z0(3:5), w0(3:5));
H3=mobius(z0(5:7), w0(5:7));
H4=mobius(z0([7;8;1]), w0([7;8;1]));

wb=[b*exp(1i*t); a*exp(1i*t)];
wb(1,g1)=evalmobius(H1, wb(1,g1));
wb(1,g2)=evalmobius(H2, wb(1,g2));
wb(1,g3)=evalmobius(H3, wb(1,g3));
wb(1,g4)=evalmobius(H4, wb(1,g4));

RHS=-A1(:,[1 end])*wb;

% Solve Laplace for the conformal mapping
w=zeros(N); w([1 end],:)=wb;
w(2:end-1,:)=sylvester(A1(2:end-1,2:end-1), A2', RHS(2:end-1,:));

% Compute jacobian determinant
J=abs(diag(r)*D*w).^2;

% Solve Poisson 
bc=zeros(2, N);
F=-ones(N);
RHS=J.*F-A1(:,[1 end])*bc;
psi=zeros(N); 
psi([1 end],:)=bc;
psi(2:end-1,:)=sylvester(A1(2:end-1,2:end-1), A2', RHS(2:end-1,:));

u=real(w);
v=imag(w);
figure(1);
surf(u(:,[1:end,1]), v(:,[1:end,1]), psi(:,[1:end,1]));
colormap(jet(256)); view(2); axis square manual; 
%shading interp; 
end