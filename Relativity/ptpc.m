function [P1,P2] = ptpc(varargin)
% Partial Trace Preconditioner
A1=varargin{1};
format=varargin{end};

switch(format)
    case 'd'
        [m,n]=size(A1);
        P1=diag(sum(A1,2));
        P2=diag(sum(A1,1));
        P1=P1/n-trace(P1)/(2*m*n)*eye(m);
        P2=P2/m-trace(P2)/(2*m*n)*eye(n);
    
    case 'dx'
        [m,n]=size(A1);
        A2=varargin{2};
        P1=bsxfun(@times, sum(A1,2), A2);
        P2=diag(A1'*diag(A2));
        P1=P1/n-trace(P1)/(2*m*n)*eye(m);
        P2=P2/m-trace(P2)/(2*m*n)*eye(n);
        
    case 'dy'
        [m,n]=size(A1);
        A2=varargin{2};
        P1=diag(A1*diag(A2));
        P2=bsxfun(@times, sum(A1,1)', A2);
        P1=P1/n-trace(P1)/(2*m*n)*eye(m);
        P2=P2/m-trace(P2)/(2*m*n)*eye(n);
        
    case 'xdx'
        A2=varargin{2};
        A3=varargin{3};
        [m,n]=size(A2);
        P1=A1*diag(sum(A2,2))*A3;
        P2=zeros(n);
        P1=P1/n;
        
    case 'ydy'
        A2=varargin{2};
        A3=varargin{3};
        [m,n]=size(A2);
        P2=A1*diag(sum(A2,1))*A3;
        P1=zeros(m);
        P2=P2/m;
        
    case 'xdy'
        A2=varargin{2};
        A3=varargin{3};
        [m,n]=size(A2);
        P1=A1*diag(A2*diag(A3));
        P2=diag(diag(A1)'*A2)*A3;
        P1=P1/n-trace(P1)/(2*m*n)*eye(m);
        P2=P2/m-trace(P2)/(2*m*n)*eye(n);
        
    case 'ydx'
        A2=varargin{2};
        A3=varargin{3};
        [m,n]=size(A2);
        P1=diag(diag(A1)'*A2')*A3;
        P2=A1*diag(A2'*diag(A3));
        P1=P1/n-trace(P1)/(2*m*n)*eye(m);
        P2=P2/m-trace(P2)/(2*m*n)*eye(n);
        
end


end

