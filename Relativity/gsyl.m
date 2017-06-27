function [X] = gsyl(A1,B1,A2,B2,C)
% Generalized Sylvester equation A1*X*B1'+A2*X*B2'=C
[S1,L1]=eig(A1,A2,'vector');
[S2,L2]=eig(B2,B1,'vector');
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2;
W1=A2*S1;
W2=B1*S2;
X=S1*((W1\C/W2.')./LL)*S2.';
end