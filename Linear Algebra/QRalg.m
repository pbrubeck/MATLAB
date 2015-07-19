function L = QRalg(A)
% Returns the eigenvalues of A.
for i=1:20
    Q=GramSchmidt(A);
    A=Q'*A*Q;
end
L=diag(A);
end