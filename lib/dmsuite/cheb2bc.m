    function [xt,D2t,D1t,phip,phim]=cheb2bc(N,g);

% Program for computing first and second derivative matrices and
% and boundary condition functions for 2 point boundary conditions
%
%  a_1 u(1)  + b_1 u'(1)  = c_1
%  a_N u(-1) + b_N u'(-1) = c_N
%
%
% INPUT 
% N        =  number of Chebyshev points in [-1,1]
% g        =  boundary condition matrix = [a_1 b_1 c_1; a_N b_N c_N]
% 
% OUTPUT  
% xt       =  Chebyshev points corresponding to rows and columns
%             of D1t and D2t
% D1t      =  1st derivative matrix incorporating bc
% D2t      =  2nd derivative matrix incorporating bc
% phip     =  1st and 2nd derivative of bc function at x=1
%             (array with 2 columns)
% phim     =  1st and 2nd derivative of bc function at x=-1 
%             (array with 2 columns)

% S.C. Reddy, J.A.C. Weideman  1998


% Get differentiation matrices

   [x,DM]=chebdif(N,2);
   D0=eye(N,N);
   D1=DM(:,:,1);
   D2=DM(:,:,2);

% extract boundary condition coefficients

   a1=g(1,1); b1=g(1,2); c1=g(1,3);
   aN=g(2,1); bN=g(2,2); cN=g(2,3);

% Case 0: Invalid boundary condition information

if ((a1==0 & b1==0) | (aN==0 & bN==0)),

   fprintf('Invalid boundary condition information (no output) \n');


elseif (b1==0 & bN==0)          % Dirichlet/Dirichlet 

   J=2:N-1;
   K=(2:N-1)';
   D1t=D1(J,K);
   D2t=D2(J,K);
   phip=c1*[D1(K,1) D2(K,1)]/a1;          % phi_+
   phim=cN*[D1(K,N) D2(K,N)]/aN;          % phi_- 
   xt=x(K);                               % node vector 

elseif (b1~=0 & bN==0),         % Dirichlet x=-1, Robin x=1

   J=2:N-1; 
   K=(1:N-1)';
   xjrow=2*sin((J-1)*pi/2/(N-1)).^2;      % 1-x_j, using trig identity
   xkcol=2*sin((K-1)*pi/2/(N-1)).^2;      % 1-x_k, using trig identity
   oner=ones(size(xkcol));                % column of ones

   fac0 = oner*(1./xjrow);                % matrix -1/(1-x_j)
   fac1 = xkcol*(1./xjrow);               % matrix (1-x_k)/(1-x_j)
   D1t = fac1.*D1(K,J)-fac0.*D0(K,J);
   D2t = fac1.*D2(K,J)-2*fac0.*D1(K,J); 

   cfac = D1(1,1)+a1/b1;                  % compute phi'_1, phi''_1
   fcol1 = -cfac*D0(K,1)+(1+cfac*xkcol).*D1(K,1);
   fcol2 = -2*cfac*D1(K,1)+(1+cfac*xkcol).*D2(K,1);
   D1t  = [fcol1 D1t];                    
   D2t  = [fcol2 D2t];                    

   phim = xkcol.*D1(K,N)/2-D0(K,N)/2;     % phi'_-, phi''_- 
   phim = cN*[phim xkcol.*D2(K,N)/2-D1(K,N)]/aN;

   
   phip= -xkcol.*D1(K,1)+D0(K,1);         % phi'_+, phi''_+ 
   phip= c1*[phip -xkcol.*D2(K,1)+2*D1(K,1)]/b1;

   xt = x(K);                             % node vector

elseif (b1==0 & bN~=0),

% Case 3: Dirichlet at x=1 and Neumann or Robin boundary x=-1.

   J=2:N-1; 
   K=(2:N)';
   xjrow=2*cos((J-1)*pi/2/(N-1)).^2;      % 1+x_j, using trig identity
   xkcol=2*cos((K-1)*pi/2/(N-1)).^2;      % 1+x_k, using trig identity
   oner=ones(size(xkcol));                % column of ones

   fac0 = oner*(1./xjrow);                % matrix 1/(1+x_j)
   fac1 = xkcol*(1./xjrow);               % matrix (1+x_k)/(1+x_j)
   D1t = fac1.*D1(K,J)+fac0.*D0(K,J);
   D2t = fac1.*D2(K,J)+2*fac0.*D1(K,J); 

   cfac = D1(N,N)+aN/bN;                  % compute phi'_N, phi''_N
   lcol1 = -cfac*D0(K,N)+(1-cfac*xkcol).*D1(K,N);
   lcol2 = -2*cfac*D1(K,N)+(1-cfac*xkcol).*D2(K,N);
   D1t  = [D1t lcol1];                    
   D2t  = [D2t lcol2];                

   phip= xkcol.*D1(K,1)/2+D0(K,1);        % compute phi'_+,phi''_+
   phip= c1*[phip xkcol.*D2(K,1)/2+D1(K,1)]/a1;

   phim= xkcol.*D1(K,N)+D0(K,N);          % compute phi'_-,phi''_-
   phim= cN*[phim xkcol.*D2(K,N)+2*D1(K,N)]/bN;

   xt = x(K);                             % node vector

elseif (b1~=0 & bN~=0),

% Case 4: Neumann or Robin boundary conditions at both endpoints. 

   J=2:N-1;
   K=(1:N)';
   xkcol0=sin((K-1)*pi/(N-1)).^2;             % 1-x_k^2 using trig identity
   xkcol1=-2*x(K);                            % -2*x_k 
   xkcol2=-2*ones(size(xkcol0));              % -2
   xjrow=1./(sin((J-1)*pi/(N-1)).^2);         % 1-x_j^2 using trig identity

   fac0=xkcol0*xjrow;
   fac1=xkcol1*xjrow;
   fac2=xkcol2*xjrow;

   D1t=fac0.*D1(K,J)+fac1.*D0(K,J);
   D2t=fac0.*D2(K,J)+2*fac1.*D1(K,J)+fac2.*D0(K,J);

   omx=sin((K-1)*pi/2/(N-1)).^2;              % (1-x_k)/2 
   opx=cos((K-1)*pi/2/(N-1)).^2;              % (1+x_k)/2

   r0=opx+(0.5+D1(1,1)+a1/b1)*xkcol0/2;       % compute phi'_1, phi''_1
   r1=0.5-(0.5+D1(1,1)+a1/b1)*x;
   r2=-0.5-D1(1,1)-a1/b1;
   rcol1=r0.*D1(K,1)+r1.*D0(K,1);
   rcol2=r0.*D2(K,1)+2*r1.*D1(K,1)+r2.*D0(K,1);

   l0=omx+(0.5-D1(N,N)-aN/bN)*xkcol0/2;       % compute phi'_N, phi''_N
   l1=-0.5+(D1(N,N)+aN/bN-0.5)*x;
   l2=D1(N,N)+aN/bN-0.5;
   lcol1=l0.*D1(K,N)+l1.*D0(K,N);
   lcol2=l0.*D2(K,N)+2*l1.*D1(K,N)+l2.*D0(K,N);

   D1t=[rcol1 D1t lcol1];
   D2t=[rcol2 D2t lcol2];

   phim1=(xkcol0.*D1(K,N)+xkcol1.*D0(K,N))/2; 
   phim2=(xkcol0.*D2(K,N)+2*xkcol1.*D1(K,N)+xkcol2.*D0(K,N))/2;
   phim=cN*[phim1 phim2]/bN;                 % compute phi'_-, phi''_-

   phip1=(-xkcol0.*D1(K,1)-xkcol1.*D0(K,1))/2;
   phip2=(-xkcol0.*D2(K,1)-2*xkcol1.*D1(K,1)-xkcol2.*D0(K,1))/2;
   phip=c1*[phip1 phip2]/b1;                 % compute phi'_+, phi''_+

   xt=x(K);                                  % node vector

end;
