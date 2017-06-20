function [TrX,TrY] = ptop(A,B,C,D,E)
% Partial Traces of general operators 
% G = kron(E,D)*diag(C(:))*kron(B,A)
TrX=D*diag(C*diag(B*E))*A;
TrY=E*diag(diag(A*D)'*C)*B;
end