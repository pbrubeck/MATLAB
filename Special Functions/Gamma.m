function y=Gamma(z)
% Computes the gamma function of a complex number.
re=real(z);
if(re<0.5)
    y=pi/(sin(pi*z)*Gamma(1-z));
else
    z=z-1;
    g=7;
    k=10;
    C=Chebyshev(2*k+1);
    p=zeros(k);
    for i=0:k-1
        p(i+1)=0;
        fact=sqrt(2/pi)*exp(g+0.5);
        for j=0:i
            p(i+1)=p(i+1)+fact*C(2*i+1,2*j+1)*(j+g+0.5)^(-j-0.5);
            fact=fact*exp(1)*(2*j+1)/2;
        end
    end
    sum=p(1)/2;
    fact=1;
    for i=2:k
        fact=fact*(z-i+2)/(z+i-1);
        sum=sum+fact*p(i);
    end
    t=z+g+0.5;
    y=sqrt(2*pi)*t^(z+0.5)*exp(-t)*sum;
end
end