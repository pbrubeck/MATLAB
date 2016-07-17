function [V,D] = QRalg(A)
% Returns the real eigenvalues and corresponding eigenvectors 
% of positive definite matrix A.
V=eye(size(A,1)); tol=eps; r=1;
while r>tol
    Q=GramSchmidt(A);
    A=Q'*A*Q;
    V=V*Q;
    r=norm(A-diag(diag(A)));
end
D=diag(A);
V=V/diag(diag(V'*V));
end