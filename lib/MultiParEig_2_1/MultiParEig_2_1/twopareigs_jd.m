function [lambda,mu,XR,YR,XL,YL,res] = twopareigs_jd(A1,B1,C1,A2,B2,C2,neig,opts)

%TWOPAREIGS_JD   Jacobi-Davidson method for two-parameter eigenvalue problem
%
% [lambda,mu,XR,YR,XL,YL,res] = TWOPAREIGS_JD(A1,B1,C1,A2,B2,C2,neig,opts)
% returns neig eigenvalues of the two-parameter eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x
% A2*y = lambda*B2*y + mu*C2*y
%
% using the Jacobi-Davidson method.
%
% Input:
%   - A1, B1, C1, A2, B2, C2 : matrices
%   - neig : number of eigenvalues (1)             
%   - opts : options (see below)
%
% Output:
%   - lambda, mu : eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - XR, YR : components of decomposable right eigenvectors (eigenvector is kron(XR(:,j),YR(:,j)), such that
%       (A1-lambda(j)*B1-mu(j)*C1)*XR(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)*YR(:,j)=0
%   - XL, YL : components of decomposable left eigenvectors (eigenvector is kron(XL(:,j),YL(:,j)), such that
%       (A1-lambda(j)*B1-mu(j)*C1)'*XL(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)'*YL(:,j)=0
%   - res : norms of residuals for all J-D steps, try plot(log10(res))
%
% Options are (default values in parenthesis):
%   - delta : convergence criterion (norm of the residual) for the outer iteration (1e-8)
%   - X0r,Y0r : initial search space (rand(n1,1) and rand(n2,1))
%   - maxsize : maximum dimension of search spaces before restart (15)
%   - minsize : dimension of search spaces after restart (5)
%   - maxsteps : maximum number of outer iteration steps (100)
%   - innersteps : number of GMRES steps for the correction equation (5)
%   - innertol : tolerance in GMRES for the correction equation (0)
%   - window : number of Ritz values with the smallest |mu| that we compute, set to 0 (default) for all Ritz values
%   - extraction : extraction method, choices are: 
%        'maxdist' : the maximum distance from the target,
%        'mindist' : the minimum distance from the target (default),
%        'minres'  : the smallest residual,
%        'maxlambda': the eigenvalue with the lambda with the largest real part
%        'minmm' : the eigenvalue with the smallest |mu|
%   - showinfo : display nothing (default), 1: all values, 2: just eigenvalues
%   - target : target for the eigenvalues ([0 0])
%   - reschange : switch to minimum residual extraction when residual norm is small - (10^(-2.5))
%   - XPr,XPl,YPl,YPr  - prohibited directions (previously computed left and right eigenvectors) - ([],[],[],[])
%   - M1, M2 : left preconditioners - inverses for A1,B1,C1 and A2,B2,C2, respectively, can also be a function_handle that returns M1*x1 or M2*x2
%   - harmonic : set to 1 to use harmonic instead of Ritz values (0) - use this for interior eigenvalue
%   - forcereal : set to 1 if you know that eigenvectors and eigenvalues are real, default is 0
%   - refine : (2) number of TRQI steps to improve the eigenpair
%
% See also: THREEPAREIGS_JD, TWOPAREIG, TWOPAREIGS, TWOPAREIGS_IRA, TWOPAREIGS_KS,
% TWOPAREIGS_SI

% References:
%
%   1) M.E. Hochstenbach, B. Plestenjak,
%      A Jacobi-Davidson type method for a right definite two-parameter eigenvalue problem,
%      SIAM J. Matrix Anal. Appl. 24 (2002), 392-410.
%
%   2) M.E. Hochstenbach, B. Plestenjak, T. Kosir,
%      A Jacobi-Davidson type method for the two-parameter eigenvalue problem,
%      SIAM J. Matrix Anal. Appl. 26 (2005), 477-497.
%
%   3) M.E. Hochstenbach, B. Plestenjak,
%      Harmonic Rayleigh-Ritz extraction for the multiparameter eigenvalue problem,
%      Electron. Trans. Numer. Anal. 29 (2008) 81-96.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<7, neig = 1;  end
if nargin<8, opts = [];  end

n1 = size(A1,1); 
n2 = size(A2,1);
 
if isfield(opts,'delta'),       delta = getfield(opts,'delta');            else delta = 1e-8;            end
if isfield(opts,'minsize'),     minsize = getfield(opts,'minsize');        else minsize = 5;             end
if isfield(opts,'maxsize'),     maxsize = getfield(opts,'maxsize');        else maxsize = 15;            end
if isfield(opts,'maxsteps'),    maxsteps = getfield(opts,'maxsteps');      else maxsteps = 100;          end
if isfield(opts,'showinfo'),    showinfo = getfield(opts,'showinfo');      else showinfo = 1;            end
if isfield(opts,'target'),      target = getfield(opts,'target');          else target = [0 0];          end
if isfield(opts,'X0r'),         X0r = getfield(opts,'Y0r');                else X0r = rand(n1,1);       end
if isfield(opts,'Y0r'),         Y0r = getfield(opts,'Y0r');                else Y0r = rand(n2,1);       end
if isfield(opts,'innersteps'),  innersteps = getfield(opts,'innersteps');  else innersteps = 5;          end
if isfield(opts,'innertol'),    innertol = getfield(opts,'innertol');      else innertol = 0;            end
if isfield(opts,'window'),      window = getfield(opts,'window');          else window = 0;              end
if isfield(opts,'extraction'),  extraction = getfield(opts,'extraction');  else extraction = 'mindist';  end
if isfield(opts,'reschange'),   reschange = getfield(opts,'reschange');    else reschange = 10^(-2.5);   end
if isfield(opts,'XPr'),         XPr = getfield(opts,'XPr');                else XPr = [];                end
if isfield(opts,'XPl'),         XPl = getfield(opts,'XPl');                else XPl = [];                end
if isfield(opts,'YPr'),         YPr = getfield(opts,'YPr');                else YPr = [];                end
if isfield(opts,'YPl'),         YPl = getfield(opts,'YPl');                else YPl = [];                end
if isfield(opts,'M1'),          M1 = getfield(opts,'M1');                  else M1 = [];                 end
if isfield(opts,'M2'),          M2 = getfield(opts,'M2');                  else M2 = [];                 end
if isfield(opts,'harmonic'),    harmonic = getfield(opts,'harmonic');      else harmonic = 0;            end
if isfield(opts,'forcereal'),   forcereal = getfield(opts,'forcereal');    else forcereal = 0;           end
if isfield(opts,'refine'),      refine = getfield(opts,'refine');          else refine = 2;              end

if ~isempty(XPl)
   B1XPl = B1'*XPl; C1XPl=C1'*XPl; 
else
   B1XPl = []; C1XPl = [];
end
if ~isempty(YPl)
   B2YPl = B2'*YPl; C2YPl = C2'*YPl; 
else
   B2YPl = []; C2YPl = [];
end

% Initial search spaces and other initialization
U1 = rgs(X0r); 
U2 = rgs(Y0r); 

initialextraction = extraction;
lambda = []; mu = [];
conv = 0; % no convergence yet
step = 1; % number of steps

AU1 = A1*U1;  BU1 = B1*U1;  CU1 = C1*U1; 
AU2 = A2*U2;  BU2 = B2*U2;  CU2 = C2*U2; 

if harmonic
    HABC1 = A1-target(1)*B1-target(2)*C1;
    HABC2 = A2-target(1)*B2-target(2)*C2;
    HABCU1 = HABC1*U1;
    HABCU2 = HABC2*U2;
    Q1 = rgs(HABCU1); 
    Q2 = rgs(HABCU2); 
end    

% selection criteria for Delta0 orthogonality
if ~isempty(XPr) 
   innerprod = (XPr'*B1XPl).*(YPr'*C2YPl)-(XPr'*C1XPl).*(YPr'*B2YPl);
   distance = min(max(abs(innerprod)));
   minInnerProd = min(distance/1.1,1);
else
   minInnerProd = Inf; 
end

res = [];

% Start of the method (loop of restart cycles)
while (conv<neig) && (step<maxsteps)
   
  while (size(U1,2)<=maxsize) && (conv<neig) && (step<maxsteps)
      
    % Step 1: We find appropriate (harmonic) Ritz value in search space
    % --------------------------------------------------------------------------------------
    % Solution of smaller problem in the search base
    
    if harmonic
        SA1 = Q1'*HABCU1;  SB1 = Q1'*BU1;    SC1 = Q1'*CU1;
        SA2 = Q2'*HABCU2;  SB2 = Q2'*BU2;    SC2 = Q2'*CU2;
    else
        SA1 = U1'*AU1;     SB1 = U1'*BU1;    SC1 = U1'*CU1;
        SA2 = U2'*AU2;     SB2 = U2'*BU2;    SC2 = U2'*CU2;
    end
    
    if window>0
        [Al,Au,AXr,AYr] = twopareigs(SA1,SB1,SC1,SA2,SB2,SC2,window);
    else
        [Al,Au,AXr,AYr] = twopareig(SA1,SB1,SC1,SA2,SB2,SC2);
    end
    
    if harmonic
        Al = Al + target(1);
        Au = Au + target(2);
    end

    noexpansion = 1;
    
    while (conv<neig) && noexpansion
        
       switch extraction
           case 'maxdist'  % Ritz value farthest to the target
               ritznorm = abs(Al-target(1)).^2+abs(Au-target(2)).^2;
               [tilda,order] = sort(-ritznorm);
           case 'minmu'  % Ritz value with the smallest |mu|
               [tilda,order] = sort(abs(Au));
           case 'mindist'  % Ritz value closest to the target 
               ritznorm = abs(Al-target(1)).^2+abs(Au-target(2)).^2;
               [tilda,order] = sort(ritznorm);
           case 'maxlambda'  % Ritz value with the largest real part
               [tilda,order] = sort(-Al);
           case 'minres'  % Ritz pair with the smallest residual               
               rn = zeros(length(Al),1);
               for kk = 1:length(Al)
                   ABCU1 = AU1-Al(kk)*BU1-Au(kk)*CU1;
                   ABCU2 = AU2-Al(kk)*BU2-Au(kk)*CU2;
                   r1r = ABCU1*AXr(:,kk); % residual r1 right
                   r2r = ABCU2*AYr(:,kk); % residual r2 right
                   rn(kk,1) = norm([r1r;r2r], 'fro'); % norm of residual
               end   
               [tilda,order] = sort(rn);
           otherwise
               error('Unknown extraction option')
        end
         
        % we select the Ritz vector that is Delta0 orthogonal to computed eigenpairs
        if isempty(XPr)
            distance = zeros(length(Al),1);
        else   
            SAXr = U1*AXr;  SAYr = U2*AYr;
            innerprod = (SAXr'*B1XPl).*(SAYr'*C2YPl)-(SAXr'*C1XPl).*(SAYr'*B2YPl);
            distance = max(abs(innerprod),[],2); 
        end
        
        izb = 1;
        while distance(order(izb))>minInnerProd
            izb = izb+1;
            if izb>length(order);  break;  end
        end
        skipped = izb-1;
   
        if izb>length(order)  
            % all directions are not enough "Delta0 orthogonal", so we pick the one with the minimal distance    
            [tilda,pos] = min(distance);
            izb = 1;
            expandrandom = 1;  % safety flag that prevents taking already computed eigenvalue
        else
            pos = order(izb);
            expandrandom = 0;
        end;
        
        dddd = distance(pos); 
        la = Al(pos);  ua = Au(pos); 
        lapr(step) = la;  uapr(step) = ua;                 % best Ritz values 
        Xrapr = U1*AXr(:,pos);  Yrapr = U2*AYr(:,pos);     % best Ritz vectors
        
        if harmonic
            % TRQ improvement
            A1Xr = A1*Xrapr;  B1Xr = B1*Xrapr;  C1Xr = C1*Xrapr; 
            A2Yr = A2*Yrapr;  B2Yr = B2*Yrapr;  C2Yr = C2*Yrapr; 
            tmpM = [Xrapr'*B1Xr Xrapr'*C1Xr; Yrapr'*B2Yr Yrapr'*C2Yr];
            tmpb = [Xrapr'*A1Xr; Yrapr'*A2Yr];
            tilda = tmpM\tmpb;
            la = tilda(1);  ua = tilda(2);       
            lapr(step) = la;  uapr(step) = ua;   % corrected Ritz values 
        end

        ABC1 = A1-la*B1-ua*C1;
        ABC2 = A2-la*B2-ua*C2;
        r1r = ABC1*Xrapr;  % residual r1 right
        r2r = ABC2*Yrapr;  % residual r2 right
        tmpres(step) = norm([r1r;r2r], 'fro'); 	% norm of residual
        if showinfo==1
            disp(sprintf('%3d skip %3d/%3d, dist %4.1e, res: %4.1e  eig: (%+4.3e% +4.3ei, %+4.3e% +4.3ei)',...
                 step,skipped,length(order),dddd,tmpres(step),real(la),imag(la),real(ua),imag(ua)))
        end
        if length(res)<step
            res(step) = tmpres(step);
        end

        if (tmpres(step)<=delta) && (expandrandom==0) % Test for convergence 
            % we have a new eigenvalue
            if forcereal
                [tilda, posmax] = max(abs(Xrapr)); Xrapr = real(Xrapr/Xrapr(posmax)); Xrapr = Xrapr/norm(Xrapr);
                [tilda, posmax] = max(abs(Yrapr)); Yrapr = real(Yrapr/Yrapr(posmax)); Yrapr = Yrapr/norm(Yrapr);
            end
            if refine>0
                [la,ua,Xrapr,Yrapr] = trqi(A1,B1,C1,A2,B2,C2,Xrapr,Yrapr,refine,eps);
                ABC1 = A1-la*B1-ua*C1;
                ABC2 = A2-la*B2-ua*C2;
            end
            XPr = [XPr Xrapr]; 
            YPr = [YPr Yrapr]; 
            % left eigenvectors are computed for the selection criteria
            [tilda, Xlapr] = min_sing_vec(ABC1,1);
            [tilda, Ylapr] = min_sing_vec(ABC2,1);
            XPl = [XPl Xlapr];
            YPl = [YPl Ylapr];
            B1XPl = B1'*XPl;  C1XPl = C1'*XPl; 
            B2YPl = B2'*YPl;  C2YPl = C2'*YPl; 
            innerprod = (Xrapr'*B1XPl).*(Yrapr'*C2YPl)-(Xrapr'*C1XPl).*(Yrapr'*B2YPl);
            distance = max(abs(innerprod),[],2);
            if distance < minInnerProd 
                minInnerProd = distance/1.1;
            end
            conv = conv+1;  
            lambda = [lambda la];  mu = [mu ua];
            noexpansion = 1;  % noexpansion, we check other Ritz vectors
            extraction = initialextraction;

            if showinfo>0
                disp(sprintf('Eig (%2d): lambda: %+5.4e %+5.4ei, mu: %+5.4e %+5.4ei, step %d, selcrit: %5.5e',conv,real(la),imag(la),real(ua),imag(ua),step,minInnerProd))
            end
        else
            % no convergence in this step
            if (tmpres(step)<reschange) && (expandrandom==0) && (~strcmp(extraction,'minres'))
                if showinfo==1,  disp('Change to min. residual'),  end
                extraction='minres';
            end
            noexpansion = 0;
        end
    end % while (conv<nreig) && noexpansion
    
    % Step 2: We find new directions dx1 and dx2 for the search space 
    % --------------------------------------------------------------------------------------
    % we use orthogonal correction equation, which is 
    % preconditioned first order correction equation P2 in paper (2).                            
    c1 = (target(1)-la)*B1*Xrapr + (target(2)-ua)*C1*Xrapr;
    c2 = (target(1)-la)*B2*Yrapr + (target(2)-ua)*C2*Yrapr;
          
    if ~isempty(M1) 
         if isa(M1,'function_handle')
             tmp1 = M1(r1r);  tmp2 = M1(c1);
         else
             tmp1 = M1*r1r;   tmp2 = M1*c1;
         end
         r1rnew = -tmp1 + tmp2*(Xrapr'*tmp1)/(Xrapr'*tmp2);
    else
         tmp2 = []; 
         r1rnew = -r1r; 
    end    
    if ~isempty(M2) 
         if isa(M2,'function_handle')
             tmp3 = M2(r2r);  tmp4 = M2(c2);
         else
             tmp3 = M2*r2r;   tmp4 = M2*c2;
         end
         r2rnew = -tmp3 + tmp4*(Yrapr'*tmp3)/(Yrapr'*tmp4);
    else
         tmp4 = [];
         r2rnew = -r2r;
    end
         
    dxr1 = gmres_jd(tmp2, Xrapr, M1, c1, Xrapr, ABC1, Xrapr, Xrapr, r1rnew, innersteps, innertol);
    dxr2 = gmres_jd(tmp4, Yrapr, M2, c2, Yrapr, ABC2, Yrapr, Yrapr, r2rnew, innersteps, innertol);

    % Step 3: We expand search space in directions dx1,dx2
    % --------------------------------------------------------------------------------------
    % dx1 and dx2 are new directions, we need orthogonal projections on orth(U1) and orth(U2)
    
    k = size(U1,2)+1;
    
    newxr1 = ExpandMGS(U1,dxr1);
    newxr2 = ExpandMGS(U2,dxr2);
    U1(:,k) = newxr1;
    U2(:,k) = newxr2;
    
    % we compute new columns and rows of matrices U1'*A1*U1,U1'*B1*U1,U1'*C1*U1,U2'*A2*U2,U2'*B2*U2,U2'*C2*U2
    AU1 = [AU1 A1*newxr1];  BU1 = [BU1 B1*newxr1];  CU1 = [CU1 C1*newxr1];
    AU2 = [AU2 A2*newxr2];  BU2 = [BU2 B2*newxr2];  CU2 = [CU2 C2*newxr2];

    if harmonic
       tmpx = HABC1*newxr1;
       tmpy = HABC2*newxr2;
       HABCU1 = [HABCU1 tmpx];
       HABCU2 = [HABCU2 tmpy];
       newh1 = ExpandMGS(Q1, tmpx);
       newh2 = ExpandMGS(Q2, tmpy);
       Q1 = [Q1 newh1];
       Q2 = [Q2 newh2];
    end
    
    step = step+1;
    
  end % while

  % Step 4: We restart J-D 
  % --------------------------------------------------------------------------------------
  % we restart with minsize last approximations that are Delta0 distant enough distances
  if (conv<neig) && (step<maxsteps)
     if showinfo==1, fprintf('--restart--'), end
     
     RestXr = [];  RestYr = [];
     chosen = 0;
     % New directions are taken from last U1
     BaseXrapr = U1(:,1:k-1)*AXr; 
     BaseYrapr = U2(:,1:k-1)*AYr;
    
     k = 1;
     while (chosen<minsize) && (k<=length(order))
         pos = order(k);  
         SAXr = BaseXrapr(:,pos);
         SAYr = BaseYrapr(:,pos);
       
         if conv==0
             distance = 0;
         else
             innerprod = (SAXr'*B1XPl).*(SAYr'*C2YPl)-(SAXr'*C1XPl).*(SAYr'*B2YPl);
             distance = max(abs(innerprod),[],2); 
         end
       
         if chosen>0
             kot = min([subspace(RestXr,SAXr),subspace(RestYr,SAYr)]);
         else
             kot = 1;
         end;
               
         if (distance<minInnerProd) && (kot>1e-4)
             chosen = chosen+1;
             RestXr = [RestXr SAXr];            
             RestYr = [RestYr SAYr];            
         end
         k = k+1;
     end

     if chosen<minsize
        if showinfo==1,  fprintf('Randomly chosen %d',minsize-chosen);  end
        RestXr = [RestXr randn(n1,minsize-chosen)];
        RestYr = [RestYr randn(n2,minsize-chosen)];
     end
     if showinfo==1,  fprintf('\n'),  end
        
     U1 = rgs(RestXr); 
     U2 = rgs(RestYr);
    
     AU1 = A1*U1; BU1 = B1*U1; CU1 = C1*U1; 
     AU2 = A2*U2; BU2 = B2*U2; CU2 = C2*U2; 
     
     if harmonic
        HABCU1 = HABC1*U1;
        HABCU2 = HABC2*U2;
        Q1 = rgs(HABCU1); 
        Q2 = rgs(HABCU2); 
     end
  end  
    
end % outer loop : while (conv<nreig) && (step<maxJDsteps)

% Results
XR = XPr;  YR = YPr;                       % right eigenvectors
XL = XPl;  YL = YPl;                       % left eigenvectors
step = step-1;

lambda = lambda(:);
mu = mu(:);

% -----------------------------------------------------------------------------------------------
% Auxiliary routine ExpandMGS
% -----------------------------------------------------------------------------------------------

function y = ExpandMGS(Q,x)

% Auxiliary routine that orthogonalizes x against the orthogonal columns 
% of U. Two steps of the repeated GS are used to maintain the stability

% Bor Plestenjak
% last revision 22.08.04

c = 0.71;
k = size(Q,2);
oldnorm = norm(x);
for j = 1:k % modified GS for additional direction
   x = x - (Q(:,j)'*x)*Q(:,j);
end
newnorm = norm(x);
if newnorm < c*oldnorm
   x = x/newnorm;
   for j = 1:k % modified GS for additional direction
      x = x - (Q(:,j)'*x)*Q(:,j);
   end
   newnorm = norm(x);
end
y = x/newnorm;

% -----------------------------------------------------------------------------------------------
% Auxiliary routine rgs (repeated Gram-Schmidt orthogonalization of A)
% -----------------------------------------------------------------------------------------------

function Q = rgs(A)

c=0.71;

[m,n] = size(A);

for i = 1:n
   Q(:,i) = A(:,i);
   oldnorm = norm(Q(:,i));
   r = oldnorm;
   for j = 1:i-1
      r = Q(:,j)'*Q(:,i);
      Q(:,i) = Q(:,i)-r*Q(:,j);
   end
   newnorm = norm(Q(:,i));
   if newnorm < c*oldnorm 
       for j = 1:i-1
          r = Q(:,j)'*Q(:,i);
          Q(:,i) = Q(:,i)-r*Q(:,j);
       end
       newnorm = norm(Q(:,i));
   end
   Q(:,i) = Q(:,i)/newnorm;
end

% -----------------------------------------------------------------------------------------------
% Auxiliary routine gmres_jd
% -----------------------------------------------------------------------------------------------

function [x, relres, k] = gmres_jd(A1x, A1y, A2, A3x, A3y, A4, A5x, A5y, b, maxit, tol)

%GMRES_JD is a version of gmres that is used in correction equations 
%   in Jacobi-Davidson type methods for two-parameter eigenvalue problems 
%   (can be applied to other similar JD methods as well)
%   The matrix vector multiplicaton in each step is 
%   y = P1* A2 * P3 * A4 * P5 * x,
%   where 
%      P1 : a) P1 is projection (I-inv(A1y'*A1x)*A1x*A1y') or 
%           b) P1=I when A1x=[] 
%      A2 : a) A2=A2 or 
%           b) A2=I when A2A=[] or
%           c) A2 is given by function_handle that evaluates A2(x)
%      P3 : a) P3 is projection (I-inv(A3y'*A3x)*A3x*A3y') or 
%           b) P3=I when A3x=[] 
%      A4 : a) A4=A4 or 
%           b) A4=I when A4A=[] or
%      P5 : a) P5 is projection (I-inv(A5y'*A5x)*A5x*A5y') or 
%           b) P5=I when A5x=[] 
%
%   function [x, relres, k] = gmres_jd(A1x, A1y, A2, A3x, A3y, A4, A5x, A5y, b, maxit, tol)
%   Assumes x0 = 0
%   In:
%      b     : righthand side
%      maxit : maximum number of steps (no restarts)
%      tol   : tolerance
%   Out:
%      x     : gmres solution
%      relres: obtained residual reduction
%      k     : number of steps performed

% Last revision: 31.08.2015
% Adapted from gmres_fast and mgs code by Michiel Hochstenbach
% Bor Plestenjak

c = 0.71; % ~= 1/sqrt(2)
H = zeros(maxit+1, maxit);
V = zeros(length(b), maxit);

if (maxit < 1) || (tol >= 1)
   relres = 1;
   x = b;
   return
end 

if ~isempty(A1x),  A1d = inv(A1y'*A1x);  else A1d = 1;  end
if ~isempty(A3x),  A3d = inv(A3y'*A3x);  else A3d = 1;  end
if ~isempty(A5x),  A5d = inv(A5y'*A5x);  else A5d = 1;  end

rho0 = norm(b);
b = b/rho0;

v = b;
Gamma = 1;
k = 0;
rho = 1;
if tol == 0
	tol0 = Inf;
else
	tol0 = 1/(tol*tol); 
end

while (k < maxit) && (rho < tol0) 
   k = k+1;
   V(:,k) = v;
  
   % multiplication by A, projections, and preconditioner
   % ----------------------------------------------------
   if ~isempty(A5x),  v = v - A5x*(A5d*(A5y'*v));  end
   if ~isempty(A4),   v = A4*v;  end
   if ~isempty(A3x),  v = v - A3x*(A3d*(A3y'*v));  end
   if ~isempty(A2)
      if isa(A2,'function_handle')
          v = A2(v);
      else
          v = A2*v;
      end
   end
   if ~isempty(A1x),  v = v - A1x*(A1d*(A1y'*v));  end
   % ----------------------------------------------------
 
   % Arnoldi step 
 
   H(k+1,k) = 1;
   l = 1;
   norm_input = norm(v);
   norm_output = 0;
   while (l <= 2) && (norm_output < c * norm_input)
      for j = 1:k
         inpr = V(:,j)'*v;
         v = v - inpr*V(:,j);    
         H(j,k) = H(j,k) + H(k+1,k) * inpr;
      end
	  norm_output = norm(v);
      v = v/norm_output;
      H(k+1,k) = H(k+1,k) * norm_output;
      l = l+1;
   end
   
   gamma = H(k+1,k);
  
   if gamma == 0 % Lucky break-down
      break
   else
      gamma = -Gamma*H(1:k,k)/gamma; 
      Gamma = [Gamma gamma];
      rho   = rho + gamma'*gamma;
   end
end
if gamma == 0   % Lucky break-down
   relres = 0;
   y = zeros(k,1);
   y(1) = 1;
   x = V(:,1:k)*(H(1:k,1:k) \ y);
else            % solve in least square sense 
   relres = 1 / sqrt(rho);
   y = zeros(k+1,1);
   y(1) = 1;
   x = V(:,1:k)*(H(1:k+1,1:k) \ y);
end

x = rho0 * x;
