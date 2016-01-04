function [f] = fpt(a,b,c,g)
% Fast Polynomial Transform, based on Potts (1998)
% Basis exchange for orthogonal polynomials that satisfy the three term
% recurrence relation: P[n](x)=(a[n]x+b[n])P[n-1](x)+c[n]P[n-2](x).
% The last input argument, g, should be the coeffients of the linear
% combination of P[n](x) and f corresponds to the coefficients of the same
% polynomial expressed as a Chebyshev series.
% The resulting Chebyshev series may be evaluated at arbitrary points via a
% non-uniform discrete cosine transform (ndct).

N=length(g)-1; t=log2(N);
gg=[g(:), zeros(N+1,1)]; gg(N+1,:)=[];
gg(N-1,1)=g(N-1)+c(N)*g(N+1);
gg(N,:)=[g(N)+b(N)*g(N+1), a(N)*g(N+1)];
ggg=zeros(size(gg));
for tau=1:t-1
   step=2^(tau+1);
   hstep=2^tau;
   ggg=reshape(ggg, [], step);
   x=cos(pi*(1:2:2*step-1)/(2*step));
   for k=0:N/step-1
       i=2*k+1; j=4*k+1;
       shift=step*k+1;
       [u11,u12,u21,u22]=calculateU(a,b,c,hstep-1,shift,x);
       ggg( i ,:) = c(shift+1)*innerP(u11, u12, gg(j+2,:), gg(j+3,:));
       ggg(i+1,:) = innerP(u21, u22, gg(j+2,:), gg(j+3,:));
       ggg( i, 1:hstep) = ggg( i ,1:hstep)+gg(j,:);
       ggg(i+1,1:hstep) = ggg(i+1,1:hstep)+gg(j+1,:);
   end
   gg=ggg;
end
%T=diag([1; ones(N-1,1)/2], 1)+diag([ones(N-1,1)/2; 1], -1);
%f=[gg(1,:)'; 0]+(a(1)*T'+b(1)*eye(N+1))*[gg(2,:)'; 0];
mid=gg(1,2:end)+b(1)*gg(2,2:end)+a(1)/2*(gg(2,1:end-1)+[gg(2,3:end), 0]);
f=[gg(1,1)+b(1)*gg(2,1)+a(1)/2*gg(2,2); mid(:); a(1)/2*gg(2,end)];
end

function [u11, u12, u21, u22] = calculateU(a,b,c,n,j,x)
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

function [P] = innerP(u1, u2, a1, a2)
% Computes polynomial vector dot product [u1; u2]'*[a1; a2]
aa1=zeros(size(u1)); aa1(1:length(a1))=a1; aa1(1)=aa1(1)*sqrt(2);
aa2=zeros(size(u2)); aa2(1:length(a2))=a2; aa2(1)=aa2(1)*sqrt(2);
P=dct(u1.*idct(aa1)+u2.*idct(aa2)); P(1)=P(1)/sqrt(2);
end