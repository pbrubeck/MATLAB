function [W, Z] = trideigs(D, E)
% Calculates the eigenvalues and eigenvectors of a symmetric tridiagonal matrix.

N=numel(D);
IL=1; IU=N;
M=N; LDZ=N;
LWORK=20*N; LIWORK=10*N; INFO=0;
VL=0.0; VU=1.0; ABSTOL=0.0;

ISUPPZ=zeros(2*M, 1);
IWORK=zeros(LIWORK, 1);
W=zeros(N, 1);
Z=zeros(LDZ, M);
WORK=zeros(1, LWORK);

out=lapack('dstevr','V','A',N,D,E,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,ISUPPZ, ...
    WORK,LWORK,IWORK,LIWORK,INFO);
W=out{12};
Z=out{13};
end