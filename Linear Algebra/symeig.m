function [varargout] = symeig(A,B)
if(nargout<2)
    JOBZ='N';
else
    JOBZ='V';
end
UPLO='L';
N=size(A,1);
LWORK=max(20,N)*N;
L=zeros(N, 1);
WORK=zeros(LWORK, 1);
INFO=0;
if(nargin==1)
    out=lapack('dsyev',JOBZ,UPLO,N,A,N,L,WORK,LWORK,INFO);
    INFO=out{9};
    if(JOBZ=='V')
        varargout{1}=out{4}; % V
        varargout{2}=out{6}; % L
    else
        varargout{1}=out{6}; % L
    end
else
    ITYPE=1;
    out=lapack('dsygv',ITYPE,JOBZ,UPLO,N,A,N,B,N,L,WORK,LWORK,INFO);
    INFO=out{12};
    if(JOBZ=='V')
        varargout{1}=out{5}; % V
        varargout{2}=out{9}; % L
    else
        varargout{1}=out{9}; % L
    end
end
if(INFO)
    warning('dsygv::INFO = %d',INFO);
end
end