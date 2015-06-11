function rts = polyRoots(P)
% Returns all the complex roots of a given polynomial 
% using Laguerre's Method.
P=polyClean(P);
n=length(P)-1;
P=P/P(n+1);
H=P;
rts=NaN(1,n);
it=0;
if(n>1)
    Q=zeros(1, n+1);
    xi=max(abs(P));
    if(all(imag(P)==0) && ~(P(n)==0 && P(n-1)==0))
        s=2*(P(n)>=0)-1;
        xi=(-P(n)+s*(n-1)*sqrt(P(n)^2-2*n/(n-1)*P(n-1)))/n;
    end
    while(n>1)
        [Q, a, dP]=LaguerreStep(P, Q, n, xi);
        xi=xi-a;
        if(abs(a)<eps)
            % Root found, deflate polynomial
            rts(n)=xi;
            P=Q(2:n);
            n=n-1;
            % Choose next target
            if(abs(dP)>eps)
                xi=-xi;
            end  
        end
        it=it+1;
    end
end
rts(1)=-P(1);
rts=polish(H, rts);
end


function [Q, a, dP] = LaguerreStep(P, Q, n, xi)
Q(n+1)=1; dP=1; r=1;
for j=n:-1:3
	Q(j)=Q(j+1)*xi+P(j);
	dP=dP*xi+Q(j);
	r=r*xi+dP;
end
Q(2)=Q(3)*xi+P(2);
Q(1)=Q(2)*xi+P(1);
dP=dP*xi+Q(2);
if(abs(Q(1))>eps)
    G=dP/Q(1);
    H=G*G-2*r/Q(1);
    s=2*(G>=0)-1;
    a=n/(G+s*sqrt((n-1)*(n*H-G*G)));
else
    a=0;
end
end


function rts = polish(P, rts)
n=length(rts);
Q=zeros(1,n+1);
for i=1:n
    a=1;
    xi=rts(i);
    while(abs(a)>eps)
        [Q, a, ~]=LaguerreStep(P, Q, n, xi);
        xi=xi-a;
    end
    rts(i)=xi;
end
end
