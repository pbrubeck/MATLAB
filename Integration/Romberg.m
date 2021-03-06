function J = Romberg(f, a, b)
% Iterative trapezoidal method.
n=1;
i=1;
h=b-a;
J=h*(f(a)+f(b))/2;
while(true)
    n=2*n;
    h=h/2;
    J=J/2+h*sum(f(a+(1:2:n-1)*h));
    m=4;
    R=zeros(1,i);
    R(1)=J;
    for j=2:i
        R(j)=(m*R(j-1)-R0(j-1))/(m-1);
        if(abs(R(j)-R(j-1))<eps)
            J=R(j);
            return;
        end
        m=m*4;
    end
    R0=R;
    i=i+1;
end
end