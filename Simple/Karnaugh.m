function [B] = Karnaugh(A)
%KARNAUGH Summary of this function goes here
%   Detailed explanation goes here
P=[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
n=sqrt(numel(A));
B=P*reshape(A,[n n])*P;
end

