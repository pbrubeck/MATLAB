function L = Laguerre(n)
% Returns the coeficient matrix of the first n Laguerre polynomials.
L=zeros(n,n);
L(1,1)=1;
L(2,1)=1;
L(2,2)=-1;
for i=1:n-2
    L(i+2,1)=1;
    for j=2:i+2
        L(i+2,j)=((2*i+1)*L(i+1,j)-L(i+1,j-1)-i*L(i,j))/(i+1);
    end
end
end