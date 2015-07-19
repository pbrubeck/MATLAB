function L = Laguerre(n,a)
% Returns the coeficient matrix of the first n Laguerre polynomials.
if(nargin<2)
  a=0;
end
L=zeros(n,n);
L(1,1)=1;
L(2,1)=1+a;
L(2,2)=-1;
for i=1:n-2
    L(i+2,1)=((2*i+1+a)*L(i+1,1)-(i+a)*L(i,1))/(i+1);
    for j=2:i+2
        L(i+2,j)=((2*i+1+a)*L(i+1,j)-L(i+1,j-1)-(i+a)*L(i,j))/(i+1);
    end
end
end
