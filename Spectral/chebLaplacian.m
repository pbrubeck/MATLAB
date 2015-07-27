function [ L ] = chebLaplacian( N )
%CHEBLAPLACIAN Summary of this function goes here
%   Detailed explanation goes here
[D, x]=chebD(N);
D2=D^2; 
D2=sparse(D2(2:N, 2:N));
I=speye(N-1);
L=kron(I,D2)+kron(D2,I);
end