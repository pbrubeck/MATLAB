function [f] = fpt(a,b,c,g)
% Fast Polynomial Transform, based on Potts (1998).
% Basis exchange for orthogonal polynomials that satisfy the three term
% recurrence relation: P[n](x)=(a[n]x+b[n])P[n-1](x)+c[n]P[n-2](x).
% The last input argument, g, should be the coeffients of the linear
% combination of P[n](x) and f corresponds to the coefficients of the same
% polynomial expressed as a Chebyshev series.
% The resulting Chebyshev series may be evaluated at arbitrary points via a
% non-uniform discrete cosine transform (ndct).
%
% Input examples:
% LegendreP k=0:2^N, a=(2*k+1)./(k+1), b=0*k, c=-k./(k+1)
% LaguerreL k=0:2^N, a=-1./(k+1), b=(2*k+1)./(k+1), c=-k./(k+1)

global fileID;
global create;
N=length(g)-1;
t=log2(N);

filename=sprintf('fptU%d.bin',N);
create=1-exist(filename,'file')/2;
if(create)
    fileID=fopen(filename,'w');
else
    fileID=fopen(filename,'r');
end

gg=[g(:), zeros(N+1,1)];
gg(N+1,:)=[];
gg(N-1,1)=g(N-1)+c(N)*g(N+1);
gg(N,:)=[g(N)+b(N)*g(N+1), a(N)*g(N+1)];
ggg=zeros(N,2);

for tau=1:t-1
   m=2^tau;
   ggg=reshape(ggg, [], 2*m);
   x=cos(pi*(1:2:4*m-1)/(4*m));
   for i=1:2:N/m-1
       j=2*i-1;
       h=m*(i-1)+1;
       [u11,u12,u21,u22]=calculateU(a,b,c,m-1,h,x);
       ggg(i  , : ) = c(h+1)*innerP(u11, u12, gg(j+2,:), gg(j+3,:));
       ggg(i+1, : ) =        innerP(u21, u22, gg(j+2,:), gg(j+3,:));
       ggg(i  ,1:m) = ggg(i  ,1:m)+gg(j  ,:);
       ggg(i+1,1:m) = ggg(i+1,1:m)+gg(j+1,:);
   end
   gg=ggg;
end

% f=gg(1,:)+(b(1)*I+a(1)*T')*gg(2,:)
f=zeros(size(g)); 
f(1)=gg(1,1)+b(1)*gg(2,1)+a(1)*gg(2,2)/2;
f(2)=gg(1,2)+b(1)*gg(2,2)+a(1)*(gg(2,1)+gg(2,3)/2);
f(3:N-1)=gg(1,3:end-1)+b(1)*gg(2,3:end-1)+a(1)/2*(gg(2,2:end-2)+gg(2,4:end));
f(end-1)=gg(1,end)+b(1)*gg(2,end)+a(1)/2*gg(2,end-1);
f(end)=a(1)*gg(2,end)/2;
fclose(fileID);
end

function [u11, u12, u21, u22] = calculateU(a,b,c,n,j,x)
% Calculates the 4x4 matrix with polynomial entries
global fileID;
global create;
if(create)
    [u12, u11]=evalP(a,b,c,n,j+1,x);
    [u22, u21]=evalP(a,b,c,n+1,j,x);
    fwrite(fileID, [u11; u12; u21; u22], 'double');
else
    U=fread(fileID, [4, length(x)], 'double');
    u11=U(1,:); u12=U(2,:); u21=U(3,:); u22=U(4,:);
end
end

function [P, PP] = evalP(a,b,c,n,j,x)
% Returns P[n](x,j) and P[n-1](x,j)
PP=zeros(size(x)); P=ones(size(x));
for i=1:n
    temp=P;
    P=(a(i+j)*x+b(i+j)).*P+c(i+j)*PP;
    PP=temp;
end
end

function [P] = innerP(u1, u2, a1, a2)
% Computes polynomial vector dot product [u1; u2]'*[a1; a2]
a1(1)=a1(1)*sqrt(2); 
a2(1)=a2(1)*sqrt(2);
P=dct(u1.*idct(a1,length(u1))+u2.*idct(a2,length(u2)));
P(1)=P(1)/sqrt(2);
end