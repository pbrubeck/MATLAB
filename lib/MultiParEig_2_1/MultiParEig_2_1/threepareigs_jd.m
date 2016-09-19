function [lambda,mu,eta,XR,YR,ZR,XL,YL,ZL,res] = threepareigs_jd(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)

%THREEPAREIGS_JD   Jacobi-Davidson method for a three-parameter eigenvalue problem
%
% [lambda,mu,eta,XR,YR,ZR,XL,YL,ZL,res] = THREEPAREIGS_JD(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)
% returns neig eigenvalues of the three-parameter eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x + eta*D1*x
% A2*y = lambda*B2*y + mu*C2*y + eta*D2*y
% A3*z = lambda*B3*z + mu*C3*z + eta*D3*z
%
% using the Jacobi-Davidson method.
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices
%   - neig : number of eigenvalues (1)             
%   - opts : options (see below)
%
% Output:
%   - lambda, mu, eta : eigenvalue parts (eigenvalues are (lambda(j),mu(j),eta(j))
%   - XR, YR, ZR : components of decomposable right eigenvectors 
%     (eigenvector is kron(XR(:,j),kron(YR(:,j),ZR(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)*XR(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)*YR(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)*ZR(:,j)=0
%   - XL, YL, ZL : components of decomposable left eigenvectors 
%     (eigenvector is kron(XL(:,j),kron(YL(:,j),ZL(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)'*XL(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)'*YL(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)'*ZL(:,j)=0
%   - res : norms of residuals for all J-D steps, try plot(log10(res))
%
% Options are (default values in parenthesis):
%   - delta : convergence criterion (norm of the residual) for the outer iteration (1e-8)
%   - X0r,Y0r,Z0r : initial search space (rand(n1,1), rand(n2,1), rand(n3,1))
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
%        'mineta' : the eigenvalue with the smallest |eta|
%   - showinfo : display nothing (default), 1: all values, 2: just eigenvalues
%   - target : target for the eigenvalues ([0 0 0])
%   - reschange : switch to minimum residual extraction when residual norm is small - (10^(-2.5))
%   - XPr,XPl,YPl,YPr,ZPr,ZPl  - prohibited directions (previously computed left and right eigenvectors) - ([],[],[],[],[],[])
%   - M1, M2, M3 : left preconditioners -  can also be a function_handle that returns M1*x1, M2*x2, or M3*x3
%   - harmonic : set to 1 to use harmonic instead of Ritz values (0) - use this for interior eigenvalue
%   - forcereal : set to 1 if you know that eigenvectors and eigenvalues are real, default is 0
%   - refine : (2) number of TRQI steps to improve the eigenpair
%
% See also: THREEPAREIG, THREEPAREIGS, THREEPAREIGS_SI, TWOPAREIGS_JD

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

if nargin<13, neig = 1;  end
if nargin<14, opts = [];  end

n1 = size(A1,1); 
n2 = size(A2,1);
n3 = size(A3,1);
 
if isfield(opts,'delta'),       delta = getfield(opts,'delta');            else delta = 1e-8;            end
if isfield(opts,'minsize'),     minsize = getfield(opts,'minsize');        else minsize = 5;             end
if isfield(opts,'maxsize'),     maxsize = getfield(opts,'maxsize');        else maxsize = 15;            end
if isfield(opts,'maxsteps'),    maxsteps = getfield(opts,'maxsteps');      else maxsteps = 100;          end
if isfield(opts,'showinfo'),    showinfo = getfield(opts,'showinfo');      else showinfo = 1;            end
if isfield(opts,'target'),      target = getfield(opts,'target');          else target = [0 0 0];        end
if isfield(opts,'X0r'),         X0r = getfield(opts,'X0r');                else X0r = rand(n1,1);        end
if isfield(opts,'Y0r'),         Y0r = getfield(opts,'Y0r');                else Y0r = rand(n2,1);        end
if isfield(opts,'Z0r'),         Z0r = getfield(opts,'Z0r');                else Z0r = rand(n3,1);        end
if isfield(opts,'innersteps'),  innersteps = getfield(opts,'innersteps');  else innersteps = 5;          end
if isfield(opts,'innertol'),    innertol = getfield(opts,'innertol');      else innertol = 0;            end
if isfield(opts,'window'),      window = getfield(opts,'window');          else window = 0;              end
if isfield(opts,'extraction'),  extraction = getfield(opts,'extraction');  else extraction = 'mindist';  end
if isfield(opts,'reschange'),   reschange = getfield(opts,'reschange');    else reschange = 10^(-2.5);   end
if isfield(opts,'XPr'),         XPr = getfield(opts,'XPr');                else XPr = [];                end
if isfield(opts,'XPl'),         XPl = getfield(opts,'XPl');                else XPl = [];                end
if isfield(opts,'YPr'),         YPr = getfield(opts,'YPr');                else YPr = [];                end
if isfield(opts,'YPl'),         YPl = getfield(opts,'YPl');                else YPl = [];                end
if isfield(opts,'ZPr'),         ZPr = getfield(opts,'ZPr');                else ZPr = [];                end
if isfield(opts,'ZPl'),         ZPl = getfield(opts,'ZPl');                else ZPl = [];                end
if isfield(opts,'M1'),          M1 = getfield(opts,'M1');                  else M1 = [];                 end
if isfield(opts,'M2'),          M2 = getfield(opts,'M2');                  else M2 = [];                 end
if isfield(opts,'M3'),          M3 = getfield(opts,'M3');                  else M3 = [];                 end
if isfield(opts,'harmonic'),    harmonic = getfield(opts,'harmonic');      else harmonic = 0;            end
if isfield(opts,'forcereal'),   forcereal = getfield(opts,'forcereal');    else forcereal = 0;           end
if isfield(opts,'refine'),      refine = getfield(opts,'refine');          else refine = 2;              end

opts3.usesparse = 0;

if ~isempty(XPl)
   B1XPl = B1'*XPl; C1XPl = C1'*XPl; D1XPl = D1'*XPl;
   B2YPl = B2'*YPl; C2YPl = C2'*YPl; D2YPl = D2'*YPl;
   B3YPl = B3'*ZPl; C3ZPl = C3'*ZPl; D3ZPl = D3'*ZPl;
else
   B1XPl = []; C1XPl = []; D1XPl = [];
   B2YPl = []; C2YPl = []; D2YPl = [];
   B3ZPl = []; C3ZPl = []; D3ZPl = [];
end
% Initial search spaces and other initialization
U1 = rgs(X0r); 
U2 = rgs(Y0r); 
U3 = rgs(Z0r); 

initialextraction = extraction;
lambda = []; mu = []; eta = [];
conv = 0; % no convergence yet
step = 1; % number of steps

AU1 = A1*U1;  BU1 = B1*U1;  CU1 = C1*U1;  DU1 = D1*U1;
AU2 = A2*U2;  BU2 = B2*U2;  CU2 = C2*U2;  DU2 = D2*U2;
AU3 = A3*U3;  BU3 = B3*U3;  CU3 = C3*U3;  DU3 = D3*U3;

if harmonic
    HABC1 = A1-target(1)*B1-target(2)*C1-target(3)*D1;
    HABC2 = A2-target(1)*B2-target(2)*C2-target(3)*D2;
    HABC3 = A3-target(1)*B3-target(2)*C3-target(3)*D3;
    HABCU1 = HABC1*U1;
    HABCU2 = HABC2*U2;
    HABCU3 = HABC3*U3;
    Q1 = rgs(HABCU1); 
    Q2 = rgs(HABCU2); 
    Q3 = rgs(HABCU3); 
end    

% selection criteria for Delta0 orthogonality
if ~isempty(XPr) 
   innerprod = (XPr'*B1XPl).*(YPr'*C2YPl).*(ZPr'*D3ZPl) + (XPr'*C1XPl).*(YPr'*D2YPl).*(ZPr'*B3ZPl) ...
             + (XPr'*D1XPl).*(YPr'*B2YPl).*(ZPr'*C3ZPl) - (XPr'*B1XPl).*(YPr'*D2YPl).*(ZPr'*C3ZPl) ...
             - (XPr'*C1XPl).*(YPr'*B2YPl).*(ZPr'*D3ZPl) - (XPr'*D1XPl).*(YPr'*C2YPl).*(ZPr'*B3ZPl);
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
        SA1 = Q1'*HABCU1;  SB1 = Q1'*BU1;  SC1 = Q1'*CU1;  SD1 = Q1'*DU1;
        SA2 = Q2'*HABCU2;  SB2 = Q2'*BU2;  SC2 = Q2'*CU2;  SD2 = Q2'*DU2;
        SA3 = Q3'*HABCU3;  SB3 = Q3'*BU3;  SC3 = Q3'*CU3;  SD3 = Q3'*DU3;
    else
        SA1 = U1'*AU1;  SB1 = U1'*BU1;  SC1 = U1'*CU1;  SD1 = U1'*DU1;
        SA2 = U2'*AU2;  SB2 = U2'*BU2;  SC2 = U2'*CU2;  SD2 = U2'*DU2;
        SA3 = U3'*AU3;  SB3 = U3'*BU3;  SC3 = U3'*CU3;  SD3 = U3'*DU3;
    end
    
    if window>0
        [Al,Au,Ae,AXr,AYr,AZr] = threepareigs(SA1,SB1,SC1,SD1,SA2,SB2,SC2,SD2,SA3,SB3,SC3,SD3,window,opts3);
    else
        [Al,Au,Ae,AXr,AYr,AZr] = threepareig(SA1,SB1,SC1,SD1,SA2,SB2,SC2,SD2,SA3,SB3,SC3,SD3);
    end
    
    if harmonic
        Al = Al + target(1);
        Au = Au + target(2);
        Ae = Ae + target(3);
    end

    noexpansion = 1;
    
    while (conv<neig) && noexpansion
        
       switch extraction
           case 'maxdist'  % Ritz value farthest to the target
               ritznorm = abs(Al-target(1)).^2 + abs(Au-target(2)).^2 + abs(Ae-target(3)).^2;
               [tilda,order] = sort(-ritznorm);
           case 'mineta'  % Ritz value with the smallest |eta|
               [tilda,order] = sort(abs(Ae));
           case 'mindist'  % Ritz value closest to the target 
               ritznorm = abs(Al-target(1)).^2 + abs(Au-target(2)).^2 + abs(Ae-target(3)).^2;
               [tilda,order] = sort(ritznorm);
           case 'maxlambda'  % Ritz value with the largest real part
               [tilda,order] = sort(-Al);
           case 'minres'  % Ritz pair with the smallest residual               
               rn = zeros(length(Al),1);
               for kk = 1:length(Al)
                   ABCU1 = AU1-Al(kk)*BU1-Au(kk)*CU1-Ae(kk)*DU1;
                   ABCU2 = AU2-Al(kk)*BU2-Au(kk)*CU2-Ae(kk)*DU2;
                   ABCU3 = AU3-Al(kk)*BU3-Au(kk)*CU3-Ae(kk)*DU3;
                   r1r = ABCU1*AXr(:,kk); % residual r1 right
                   r2r = ABCU2*AYr(:,kk); % residual r2 right
                   r3r = ABCU3*AZr(:,kk); % residual r3 right
                   rn(kk,1) = norm([r1r;r2r;r3r], 'fro'); % norm of residual
               end   
               [tilda,order] = sort(rn);
           otherwise
               error('Unknown extraction option')
        end
         
        % we select the Ritz vector that is Delta0 orthogonal to computed eigenpairs
        if isempty(XPr)
            distance = zeros(length(Al),1);
        else   
            SAXr = U1*AXr;  SAYr = U2*AYr;  SAZr = U3*AZr;
            innerprod = (SAXr'*B1XPl).*(SAYr'*C2YPl).*(SAZr'*D3ZPl) + (SAXr'*C1XPl).*(SAYr'*D2YPl).*(SAZr'*B3ZPl) ...
                      + (SAXr'*D1XPl).*(SAYr'*B2YPl).*(SAZr'*C3ZPl) - (SAXr'*B1XPl).*(SAYr'*D2YPl).*(SAZr'*C3ZPl) ...
                      - (SAXr'*C1XPl).*(SAYr'*B2YPl).*(SAZr'*D3ZPl) - (SAXr'*D1XPl).*(SAYr'*C2YPl).*(SAZr'*B3ZPl);
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
        la = Al(pos);  ua = Au(pos);  ea = Ae(pos);
        lapr(step) = la;  uapr(step) = ua;  eapr(step) = ea;  % best Ritz values 
        Xrapr = U1*AXr(:,pos);  Yrapr = U2*AYr(:,pos); Zrapr = U3*AZr(:,pos);  % best Ritz vectors
        
        if harmonic
            % TRQ improvement
            A1Xr = A1*Xrapr;  B1Xr = B1*Xrapr;  C1Xr = C1*Xrapr;  D1Xr = D1*Xrapr;  
            A2Yr = A2*Yrapr;  B2Yr = B2*Yrapr;  C2Yr = C2*Yrapr;  D2Yr = D2*Yrapr;  
            A3Zr = A3*Zrapr;  B3Zr = B3*Zrapr;  C3Zr = C3*Zrapr;  D3Zr = D3*Zrapr;  
            tmpM = [Xrapr'*B1Xr Xrapr'*C1Xr Xrapr'*D1Xr; Yrapr'*B2Yr Yrapr'*C2Yr Yrapr'*D2Yr; Zrapr'*B3Zr Zrapr'*C3Zr Zrapr'*D3Zr];
            tmpb = [Xrapr'*A1Xr; Yrapr'*A2Yr;  Zrapr'*A3Zr];
            tilda = tmpM\tmpb;
            la = tilda(1);  ua = tilda(2); ea = tilda(3);      
            lapr(step) = la;  uapr(step) = ua;  eapr(step) = ea; % corrected Ritz values 
        end

        ABC1 = A1-la*B1-ua*C1-ea*D1;
        ABC2 = A2-la*B2-ua*C2-ea*D2;
        ABC3 = A3-la*B3-ua*C3-ea*D3;
        r1r = ABC1*Xrapr;  % residual r1 right
        r2r = ABC2*Yrapr;  % residual r2 right
        r3r = ABC3*Zrapr;  % residual r3 right
        tmpres(step) = norm([r1r;r2r;r3r], 'fro'); 	% norm of residual
        if showinfo==1
            disp(sprintf('%3d skip %3d/%3d, dist %4.1e, res: %4.1e  eig: (%+4.3e% +4.3ei, %+4.3e% +4.3ei, %+4.3e% +4.3ei)',...
                 step,skipped,length(order),dddd,tmpres(step),real(la),imag(la),real(ua),imag(ua),real(ea),imag(ea)))
        end
        if length(res)<step
            res(step) = tmpres(step);
        end

        if (tmpres(step)<=delta) && (expandrandom==0) % Test for convergence 
            % we have a new eigenvalue
            if forcereal
                [tilda, posmax] = max(abs(Xrapr)); Xrapr = real(Xrapr/Xrapr(posmax)); Xrapr = Xrapr/norm(Xrapr);
                [tilda, posmax] = max(abs(Yrapr)); Yrapr = real(Yrapr/Yrapr(posmax)); Yrapr = Yrapr/norm(Yrapr);
                [tilda, posmax] = max(abs(Zrapr)); Zrapr = real(Zrapr/Zrapr(posmax)); Zrapr = Zrapr/norm(Zrapr);
            end
            if refine>0
                [la,ua,ea,Xrapr,Yrapr,Zrapr] = trqi_3p(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,Xrapr,Yrapr,Zrapr,refine,eps);
                ABC1 = A1-la*B1-ua*C1-ea*D1;
                ABC2 = A2-la*B2-ua*C2-ea*D2;
                ABC3 = A3-la*B3-ua*C3-ea*D3;
            end
            XPr = [XPr Xrapr]; 
            YPr = [YPr Yrapr]; 
            ZPr = [ZPr Zrapr]; 
            % left eigenvectors are computed for the selection criteria
            [tilda, Xlapr] = min_sing_vec(ABC1,1);
            [tilda, Ylapr] = min_sing_vec(ABC2,1);
            [tilda, Zlapr] = min_sing_vec(ABC3,1);
            XPl = [XPl Xlapr];
            YPl = [YPl Ylapr];
            ZPl = [ZPl Zlapr];
            B1XPl = B1'*XPl; C1XPl = C1'*XPl; D1XPl = D1'*XPl;
            B2YPl = B2'*YPl; C2YPl = C2'*YPl; D2YPl = D2'*YPl;
            B3ZPl = B3'*ZPl; C3ZPl = C3'*ZPl; D3ZPl = D3'*ZPl;
            innerprod = (Xrapr'*B1XPl).*(Yrapr'*C2YPl).*(Zrapr'*D3ZPl) + (Xrapr'*C1XPl).*(Yrapr'*D2YPl).*(Zrapr'*B3ZPl) ...
                      + (Xrapr'*D1XPl).*(Yrapr'*B2YPl).*(Zrapr'*C3ZPl) - (Xrapr'*B1XPl).*(Yrapr'*D2YPl).*(Zrapr'*C3ZPl) ...
                      - (Xrapr'*C1XPl).*(Yrapr'*B2YPl).*(Zrapr'*D3ZPl) - (Xrapr'*D1XPl).*(Yrapr'*C2YPl).*(Zrapr'*B3ZPl);
            distance = max(abs(innerprod),[],2);
            if distance < minInnerProd 
                minInnerProd = distance/1.1;
            end
            conv = conv+1;  
            lambda = [lambda la];  mu = [mu ua];  eta = [eta ea];
            noexpansion = 1;  % noexpansion, we check other Ritz vectors
            extraction = initialextraction;

            if showinfo>0
                disp(sprintf('Eig (%2d): lambda: %+5.4e %+5.4ei, mu: %+5.4e %+5.4ei, eta: %+5.4e %+5.4ei, step %d, selcrit: %5.5e',conv,real(la),imag(la),real(ua),imag(ua),real(ea),imag(ea),step,minInnerProd))
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
    c1 = (target(1)-la)*B1*Xrapr + (target(2)-ua)*C1*Xrapr + (target(3)-ea)*D1*Xrapr;
    c2 = (target(1)-la)*B2*Yrapr + (target(2)-ua)*C2*Yrapr + (target(3)-ea)*D2*Yrapr;
    c3 = (target(1)-la)*B3*Zrapr + (target(3)-ua)*C3*Zrapr + (target(3)-ea)*D3*Zrapr;
          
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
    if ~isempty(M3) 
         if isa(M3,'function_handle')
             tmp5 = M3(r3r);  tmp6 = M3(c3);
         else
             tmp5 = M3*r3r;   tmp6 = M3*c3;
         end
         r3rnew = -tmp5 + tmp6*(Zrapr'*tmp5)/(Zrapr'*tmp6);
    else
         tmp6 = [];
         r3rnew = -r3r;
    end
         
    dxr1 = gmres_jd(tmp2, Xrapr, M1, c1, Xrapr, ABC1, Xrapr, Xrapr, r1rnew, innersteps, innertol);
    dxr2 = gmres_jd(tmp4, Yrapr, M2, c2, Yrapr, ABC2, Yrapr, Yrapr, r2rnew, innersteps, innertol);
    dxr3 = gmres_jd(tmp6, Zrapr, M3, c3, Zrapr, ABC3, Zrapr, Zrapr, r3rnew, innersteps, innertol);

    % Step 3: We expand search space in directions dx1,dx2,dx3
    % --------------------------------------------------------------------------------------
    % dx1, dx2, dx3 are new directions, we need orthogonal projections on
    % orth(U1), orth(U2), and orth(U3)
    
    k = size(U1,2)+1;
    
    newxr1 = ExpandMGS(U1,dxr1);
    newxr2 = ExpandMGS(U2,dxr2);
    newxr3 = ExpandMGS(U3,dxr3);
    U1(:,k) = newxr1;
    U2(:,k) = newxr2;
    U3(:,k) = newxr3;
    
    % we compute new columns and rows of matrices
    % U1'*A1*U1,U1'*B1*U1,U1'*C1*U1,U1'*D1*U1, U2'*A2*U2,U2'*B2*U2,U2'*C2*U2, ..,
    % U3'*D3*U3
    AU1 = [AU1 A1*newxr1];  BU1 = [BU1 B1*newxr1];  CU1 = [CU1 C1*newxr1];  DU1 = [DU1 D1*newxr1];
    AU2 = [AU2 A2*newxr2];  BU2 = [BU2 B2*newxr2];  CU2 = [CU2 C2*newxr2];  DU2 = [DU2 D2*newxr2];
    AU3 = [AU3 A3*newxr3];  BU3 = [BU3 B3*newxr3];  CU3 = [CU3 C3*newxr3];  DU3 = [DU3 D3*newxr3];

    if harmonic
       tmpx = HABC1*newxr1;
       tmpy = HABC2*newxr2;
       tmpz = HABC3*newxr3;
       HABCU1 = [HABCU1 tmpx];
       HABCU2 = [HABCU2 tmpy];
       HABCU3 = [HABCU3 tmpz];
       newh1 = ExpandMGS(Q1, tmpx);
       newh2 = ExpandMGS(Q2, tmpy);
       newh3 = ExpandMGS(Q3, tmpz);
       Q1 = [Q1 newh1];
       Q2 = [Q2 newh2];
       Q3 = [Q3 newh3];
    end
    
    step = step+1;
    
  end % while

  % Step 4: We restart J-D 
  % --------------------------------------------------------------------------------------
  % we restart with minsize last approximations that are Delta0 distant enough distances
  if (conv<neig) && (step<maxsteps)
     if showinfo==1, fprintf('--restart--'), end
     
     RestXr = [];  RestYr = [];  RestZr = [];
     chosen = 0;
     % New directions are taken from last U1
     BaseXrapr = U1(:,1:k-1)*AXr; 
     BaseYrapr = U2(:,1:k-1)*AYr;
     BaseZrapr = U3(:,1:k-1)*AZr;
    
     k = 1;
     while (chosen<minsize) && (k<=length(order))
         pos = order(k);  
         SAXr = BaseXrapr(:,pos);
         SAYr = BaseYrapr(:,pos);
         SAZr = BaseZrapr(:,pos);
       
         if conv==0
             distance = 0;
         else
             innerprod = (SAXr'*B1XPl).*(SAYr'*C2YPl).*(SAZr'*D3ZPl) + (SAXr'*C1XPl).*(SAYr'*D2YPl).*(SAZr'*B3ZPl) ...
                       + (SAXr'*D1XPl).*(SAYr'*B2YPl).*(SAZr'*C3ZPl) - (SAXr'*B1XPl).*(SAYr'*D2YPl).*(SAZr'*C3ZPl) ...
                       - (SAXr'*C1XPl).*(SAYr'*B2YPl).*(SAZr'*D3ZPl) - (SAXr'*D1XPl).*(SAYr'*C2YPl).*(SAZr'*B3ZPl);
             distance = max(abs(innerprod),[],2); 
         end
       
         if chosen>0
             kot = min([subspace(RestXr,SAXr),subspace(RestYr,SAYr),subspace(RestZr,SAZr)]);
         else
             kot = 1;
         end;
               
         if (distance<minInnerProd) && (kot>1e-4)
             chosen = chosen+1;
             RestXr = [RestXr SAXr];            
             RestYr = [RestYr SAYr];            
             RestZr = [RestZr SAZr];            
         end
         k = k+1;
     end

     if chosen<minsize
        if showinfo==1,  fprintf('Randomly chosen %d',minsize-chosen);  end
        RestXr = [RestXr randn(n1,minsize-chosen)];
        RestYr = [RestYr randn(n2,minsize-chosen)];
        RestZr = [RestZr randn(n3,minsize-chosen)];
     end
     if showinfo==1,  fprintf('\n'),  end
        
     U1 = rgs(RestXr); 
     U2 = rgs(RestYr);
     U3 = rgs(RestZr);
    
     AU1 = A1*U1;  BU1 = B1*U1;  CU1 = C1*U1;  DU1 = D1*U1;
     AU2 = A2*U2;  BU2 = B2*U2;  CU2 = C2*U2;  DU2 = D2*U2;
     AU3 = A3*U3;  BU3 = B3*U3;  CU3 = C3*U3;  DU3 = D3*U3;
     
     if harmonic
        HABCU1 = HABC1*U1;
        HABCU2 = HABC2*U2;
        HABCU3 = HABC3*U3;
        Q1 = rgs(HABCU1); 
        Q2 = rgs(HABCU2); 
        Q3 = rgs(HABCU3); 
     end
  end  
    
end % outer loop : while (conv<nreig) && (step<maxJDsteps)

% Results
XR = XPr;  YR = YPr;  ZR = ZPr;            % right eigenvectors
XL = XPl;  YL = YPl;  ZL = ZPl;            % left eigenvectors
step = step-1;

lambda = lambda(:);
mu = mu(:);
eta = eta(:);

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
