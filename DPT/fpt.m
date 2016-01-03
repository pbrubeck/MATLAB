function [f] = fpt(a,b,c,g)
% Fast Polynomial Transform, based on Potts (1998)
% Basis exchange for orthogonal polynomials that satisfy the three term
% recurrence relation: P[n+1](x)=(a[n]x+b[n])P[n](x)+c[n]P[n-1](x).
% The last input argument, g, should be the coeffients of the linear
% combination of P[n](x) and f corresponds to the coefficients of the same
% polynomial expressed as a Chebyshev series.
% The resulting Chebyshev series may be evaluated at arbitrary points via a
% non-uniform discrete cosine transform (ndct).

N=length(g)-1;
gg=[g(:), zeros(N+1,1)]; gg(N+1,:)=[];
gg(N-1,1)=g(N-1)+c(N)*g(N+1);
gg(N,1)=g(N)+b(N)*g(N+1);
gg(N,2)=a(N)*g(N+1);
t=log2(N);
for tau=1:t-1
   n=2^tau-1;
   step=2^(tau+1);
   hstep=2^tau;
   ggg=zeros(N+1,step);
   for j=1:step:N+1-step
       [u11,u12,u21,u22]=calculateU(a,b,c,n,j,step);
       ggg(j,:)=c(j+1)*(fpm(gg(j+hstep,:),u11)+fpm(gg(j+hstep+1,:),u12));
       ggg(j,1:hstep)=ggg(j,1:hstep)+gg(j,:);
       ggg(j+1,:)=fpm(gg(j+hstep,:),u21)+fpm(gg(j+hstep+1,:),u22);
       ggg(j+1,1:hstep)=ggg(j+1,1:hstep)+gg(j+1,:);
   end
   gg=ggg;
end
T=diag([1; ones(N-1,1)/2], 1)+diag([ones(N-1,1)/2; 1], -1);
f=[gg(1,:)'; 0]+(a(1)*T'+b(1)*eye(N+1))*[gg(2,:)'; 0];
end

function [u11, u12, u21, u22] = calculateU(a,b,c,n,j,M)
x=cos(pi*(1:2:2*M-1)/(2*M));
[u12, u11]=evalP(a,b,c,n,j+1,x);
[u22, u21]=evalP(a,b,c,n+1,j,x);
end

function [P, PP] = evalP(a,b,c,n,j,x)
% Returns P[n](x,j) and P[n-1](x,j)
P=ones(size(x)); PP=zeros(size(x));
for i=1:n
    temp=P;
    P=(a(i+j)*x+b(i+j)).*P+c(i+j)*PP;
    PP=temp;
end
end