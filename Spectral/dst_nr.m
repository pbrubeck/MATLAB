N=16;
f0=linspace(1, N-1, N-1)';

f=[0; f0];
y=zeros(size(f));
for j=0:N-1
    Nmj=mod(N-j,N);
    y(j+1)=sin(pi*j/N)*(f(j+1)+f(Nmj+1))+0.5*(f(j+1)-f(Nmj+1));
end


yhat=1i*fft(y);
yhat=yhat(1:N/2);
yhat(1)=real(yhat(1))+0.5i*imag(yhat(1));

F=zeros(size(f));
F(1:2:end)=real(yhat);
F(2:2:end)=cumsum(imag(yhat));
disp([F, [0; dst(f0)]]);




fprintf('\n');

f=F;
y=zeros(size(f));
for j=0:N-1
    Nmj=mod(N-j,N);
    y(j+1)=sin(pi*j/N)*(f(j+1)+f(Nmj+1))+0.5*(f(j+1)-f(Nmj+1));
end


yhat=1i*fft(y);
yhat=yhat(1:N/2);
yhat(1)=real(yhat(1))+0.5i*imag(yhat(1));

F=zeros(size(f));
F(1:2:end)=real(yhat);
F(2:2:end)=cumsum(imag(yhat));
disp(2/N*[F, [0; dst(dst(f0))]]);