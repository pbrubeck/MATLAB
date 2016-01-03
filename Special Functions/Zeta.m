function y=Zeta(s)
%Computes the Riemman Zeta function by the Lanczos approximation
b=real(s)<0;
s(b)=1-s(b);
n=max(32, max(floor(real(s(:)))));
m=ones(1,n+1);
sum=1;
for i=0:n-1
    m(i+2)=4*(n-i)*(n+1)*m(i+1)/((2*i+1)*(2*i+2));
    sum=sum+m(i+2);
end
dk=sum-1;
sgn=1;
z=zeros(size(s));
for k=1:n
    z=z+sgn*dk*k.^(-s);
    dk=dk-m(k+1);
    sgn=-sgn;
end
y=z./(sum*(1-2.^(1-s)));
% Analytic continuation
w=1-s(b);
y(b)=Gamma(s(b)).*y(b).*sin(pi/2*w)/pi.*(2*pi).^w; 
end
