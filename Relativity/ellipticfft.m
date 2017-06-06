function [A2,green,ps,kd] = ellipticfft(A1,L2,B1,rd1)
% Kept degrees of freedom
m=size(A1,1);
kd1=1:m; kd1(rd1)=[];
kd=@(uu) uu(kd1,:);

% Give-back matrix
G1=-B1(:,rd1)\B1(:,kd1);
% Schur complement
S1=A1(kd1,kd1)+A1(kd1,rd1)*G1;
% Poincare-Steklov operator
N1=zeros(m,length(rd1));
N1(rd1,:)=inv(B1(:,rd1));
ps=@(b1) N1*b1;

% Eigenfunctions
V1=zeros(m,length(kd1));
[V1(kd1,:),L1]=eig(S1,'vector');
V1(rd1,:)=G1*V1(kd1,:);
[L1,LL2]=ndgrid(L1,L2);
LL=L1+LL2; LL(abs(LL)<1e-9)=inf;
W1=inv(V1(kd1,:));

A2=@(uu) ifft(bsxfun(@times, fft(uu,[],2), L2),[],2);
green=@(uu) ifft(V1*((W1*fft(uu,[],2))./LL),[],2);
end