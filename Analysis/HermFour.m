function [] = HermFour(f, n, m)
[x, w]=GaussHermite(n); ww=log(sqrt(2)*w);

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
    log_scale=round(log(abs(h_i)));
    scale=exp(-log_scale);
    [h_i_2, h_i_1]=deal(h_i_1.*scale, h_i.*scale);
    sum_log_scale=sum_log_scale+log_scale;
end

F=H'*diag(1i.^(0:m-1))*W;
G=conj(F);
[xx,yy]=meshgrid(sqrt(2)*x);
y=f(xx, yy);
u=F*y*F.';
u=G*(-u./(xx.^2+yy.^2))*G.';

figure(1);
mesh(xx, yy, real(u));
colormap(jet);
end