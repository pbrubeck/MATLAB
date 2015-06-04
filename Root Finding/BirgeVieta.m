function rts = BirgeVieta(P)
% Returns all the real roots of a given polynomial.
P=polyClean(P);
n=length(P);
P=P/P(n);
i=0;
xi=-P(1);
Q=zeros(1,n);
rts=NaN(1,n-1);
while(n>2 && i<40)
    Q(n)=q=r=1;
    for j=n-1:-1:3
        Q(j)=Q(j+1)*xi+P(j);
        q=q*xi+Q(j);
        r=r*xi+q;
    end
    Q(2)=Q(3)*xi+P(2);
    p=Q(2)*xi+P(1);
    b=abs(p)>eps;
    if(b)
        q=q*xi+Q(2);
        num=p*q;
        den=q*q-p*r;
        if(num==0 || den==0)
            xi=xi-rand();
        else
            delta=num/den;
            xi=xi-delta;
            b=abs(delta)>eps;
        end
    end
    i=i+1;
    if(~b)
        P=Q(2:n);
        n=n-1;
        rts(n)=xi;
        i=0;
    end
end
if(n==2)
    rts(1)=-P(1);
end
rts=sort(rts);
end
