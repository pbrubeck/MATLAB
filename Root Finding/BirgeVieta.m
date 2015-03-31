function roots = BirgeVieta(P)
% Returns all the real roots of a given polynomial.
n=length(P);
if(P(n)~=1)
    for i=1:n
      P(i)=P(i)/P(n);  
    end
end
xi=P(2)/P(1);
ii=0;
Q=zeros(1,n);
R=zeros(1,n);
roots=zeros(1,n-1);
while(n>1)
    Q(n)=1;
    for i=n-1:-1:1
        Q(i)=Q(i+1)*xi+P(i);
    end
    if(abs(Q(1))>1E-15)
        R(n)=1;
        for i=n-1:-1:2
            R(i)=R(i+1)*xi+Q(i);
        end
        xi=xi-Q(1)/R(2);
    else
        P=Q(2:n);
        n=n-1;
        Q=zeros(1,n);
        R=zeros(1,n);
        roots(n)=xi;
    end
    ii=ii+1;
end
end