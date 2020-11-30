function [D, x, w] = legD(N,itype)
% Legendre spectral differentiation matrix
% CG (GLL nodes) if itype==0 (default)
% DG (GL nodes) if otherwise
if(nargin<2), itype=0; end 

if(itype==0)
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
    %p=[1; ((-1).^(1:N-2)./abs(V(1,:)))'.*sqrt((1-x(2:N-1).^2)*3/(2*N*(N-1))); -(-1)^N];
else
    k=1:N-1;
    E=k./sqrt(4*k.*k-1);
    [x,V]=trideigs(zeros(1,N), E); 
    w=2*V(1,:).^2;
end
x=x(:);
w=w(:);
D=repmat(x, [1, N]);
D=D-D'+eye(N);
p=prod(D,2);
D=(p*(1./p)')./D;
D=D-diag(sum(D,2));
end