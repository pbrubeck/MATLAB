function [sem] = MathieuS(m, n, z, q)
j=mod(m,2);
k=(m+j)/2;
E=q*ones(1,n-1);
if j==0
    D=(2:2:2*n).^2;
else
    D=[1-q, (3:2:2*n-1).^2];  
end
[~,B]=trideigs(D,E);
B=B(:,k);
BB=zeros(1,2*n+1);
BB(3-j:2:end)=B;
sem=-imag(fft(BB,numel(z)));
end