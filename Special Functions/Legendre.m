function P = Legendre(n)
%Returns the coeficient matrix of the first n Legendre Polinomials
P=zeros(n,n);
P(1,1)=1;
P(2,2)=1;
for i=1:n-2
    P(i+2,1)=-i*P(i,1)/(i+1);
    for j=2:i+2
        P(i+2,j)=((2*i+1)*P(i+1,j-1)-i*P(i,j))/(i+1);
    end
end
end