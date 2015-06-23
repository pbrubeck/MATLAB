function rts = BirgeVieta(P)
% Returns all the real roots of a given polynomial.
P=polyClean(P);
n=length(P);
P=P/P(n);
i=0;
x=-P(1);
Q=zeros(1,n);
rts=NaN(1,n-1);
while(n>2 && i<40)
    Q(n)=1;q=1;r=1;
    for j=n-1:-1:3
        Q(j)=Q(j+1)*x+P(j);
        q=q*x+Q(j);
        r=r*x+q;
    end
    Q(2)=Q(3)*x+P(2);
    p=Q(2)*x+P(1);
    b=abs(p)>eps;
    if(b)
        q=q*x+Q(2);
        num=p*q;
        den=q*q-p*r;
        if(num==0 || den==0)
            x=x-rand();
        else
            delta=num/den;
            x=x-delta;
            b=abs(delta)>eps;
        end
    end
    i=i+1;
    if(~b)
        P=Q(2:n);
        n=n-1;
        rts(n)=x;
        i=0;
    end
end
if(n==2)
    rts(1)=-P(1);
end
rts=sort(rts);
end
