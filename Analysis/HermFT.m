function [] = HermFT(f, n, m)
[x, w]=GaussHermite(n,0,1); ww=log(sqrt(2)*w);
H=zeros(m, n);
H(1,:)=pi^(-1/4)*exp(-x.^2);
H(2,:)=H(1,:).*(2*x);
W=zeros(m, n);
W(1,:)=pi^(-1/4)*sqrt(2)*w;
W(2,:)=W(1,:).*(2*x);

h_i_2=pi^(-1/4)*ones(size(x));
h_i_1=pi^(-1/4)*(2*x);

sum_log_scale=zeros(size(x));
for i=1:m-2
    h_i=2/sqrt(i+1)*x.*h_i_1-sqrt(i/(i+1))*h_i_2;
    H(i+2,:)=h_i.*exp(-x.^2+sum_log_scale);
    W(i+2,:)=h_i.*exp(ww+sum_log_scale);
    log_scale=ceil(log(abs(h_i)));
    scale=exp(-log_scale);
    [h_i_2, h_i_1]=deal(h_i_1.*scale, h_i.*scale);
    sum_log_scale=sum_log_scale+log_scale;
end


F=H'*diag(1i.^(0:m-1))*W;
G=conj(F);
[kx, ky]=meshgrid(sqrt(2)*x);
Y=f(kx, ky); Z=F*Y*F.';
u=G*(-Z./(kx.^2+ky.^2))*G.';
Ex=G*(1i*kx./(kx.^2+ky.^2).*Z)*G.';
Ey=G*(1i*ky./(kx.^2+ky.^2).*Z)*G.';
figure(1);
mesh(kx, ky, real(u));
quiver(kx, ky, real(Ex), real(Ey), 'Linewidth', 2);
colormap(jet);
end