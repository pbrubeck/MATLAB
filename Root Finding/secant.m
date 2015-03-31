function b = secant(f, a, b)
% Solves f(x)=0 for x, requires two initial aproximations.
ii=0;
yb=0;
ya=f(a);
while(ii==0 || abs(yb)>1E-15)
    yb=f(b);
    fprintf('i=%d \t a=%f \t b=%f \t f(b)=%f \n', ii, a, b, yb);
    temp=b;
    b=b-yb*(b-a)/(yb-ya);
    a=temp;
    ya=yb;
    ii=ii+1;
end
end