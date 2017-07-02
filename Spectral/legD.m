function [D, x, w] = legD(N)
% Legendre spectral differentiation matrix
x=zeros(1,N); x([1,N])=[-1,1];
w=zeros(1,N); w([1,N])=2/(N*(N-1));
if(N==2)
    D=[-1,1;-1,1]/2;
    x=x(:); w=w(:);
    return;
end

k=1:N-3;
E=sqrt((k.*(k+2))./((2*k+1).*(2*k+3)));
[x(2:N-1),V]=trideigs(zeros(1,N-2), E);
w(2:N-1)=4/3*V(1,:).^2./(1-x(2:N-1).^2);
x=x(:);
w=w(:);

X=repmat(x, [1, N]);
p=[1; ((-1).^(1:N-2)./abs(V(1,:)))'.*sqrt((1-x(2:N-1).^2)*3/(2*N*(N-1))); -(-1)^N];
dX=X-X'+eye(N);
D=(p*(1./p)')./dX;
D=D-diag(sum(D,2));
end