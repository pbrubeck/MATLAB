function H = Hermite(n)
% Returns the coeficient matrix of the first n Hermite polynomials.
H=zeros(n,n);
H(1,1)=1;
H(2,2)=2;
for i=1:n-2
    H(i+2,1)=-2*i*H(i,1);
    for j=2:i+2
        H(i+2,j)=2*H(i+1,j-1)-2*i*H(i,j);
    end
end
end