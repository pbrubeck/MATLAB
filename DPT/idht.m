function [z] = idht(f, N)
% Inverse Discrete Hermite Transform
x=sqrt(pi/(2*N))*(-N:2:N-2);
H=zeros(N,N);
H(1,:)=pi^-(1/4)*exp(-x.^2/2);
H(2,:)=sqrt(2)*x.*H(1,:);
for i=1:N-2
    H(i+2,:)=sqrt(2/(i+1))*x.*H(i+1,:)-sqrt(i/(i+1))*H(i,:);
end
z=sqrt(2*pi/N)*H*f(x');
end