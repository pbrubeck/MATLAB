function [sem] = MathieuS(m, n, z, q)
j=mod(m,2);
E=q*ones(1,n-1);
if j==0
    D=(2:2:2*n).^2;
else
    D=[1-q, (3:2:2*n-1).^2];  
end
[~,B]=trideigs(D,E);
B=B(:,(m+j)/2);
sem=zeros(size(z));
for i=0:n-1
    sem=sem+B(i+1)*sin((2*i+2-j)*z);
end
end