function y=Zeta(s)
%Computes the Riemman Zeta function by the Lanczos approximation

if(real(s)<0)
    w=1-s;
    w=Gamma(w)*Zeta(w)*sin(pi/2*s);
    y=w/pi*(2*pi)^s;
else
    n=max(25, floor(real(s)));
    m=zeros(n+1);
    m(1)=1;
    sum=1;
    for i=0:n-1
        m(i+2)=4*(n-i)*(n+1)*m(i+1)/((2*i+1)*(2*i+2));
        sum=sum+m(i+2);
    end
    dk=sum-1;
    sgn=1;
    z=0;
    for k=1:n
        z=z+sgn*dk*k^(-s);
        dk=dk-m(k+1);
        sgn=-sgn;
    end
    y=z/(sum*(1-2^(1-s)));
end

end

