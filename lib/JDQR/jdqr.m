function varargout=jdqr(varargin)
%JDQR computes eigenpairs of a square matrix or operator.
% 
%  Lambda = JDQR(A)
%  returns the absolute largest eigenvalues in a K vector Lambda. Here 
%  K=min(5,N) (unless K has been specified), where N=size(A,1). JDQR(A) 
%  (without output argument) displays the K eigenvalues.
%
%  [X,Lambda] = JDQR(A) 
%  returns the eigenvectors in the N by K matrix X and the eigenvalues
%  in the K by K diagonal matrix Lambda: A*X = X*Lambda. Lambda contains
%  the Jordan structure if there are multiple eigenvalues.
%
%  [X,Lambda,Q,S] = JDQR(A) 
%  returns also a partial Schur decomposition: S is an K by K upper
%  triangular matrix and Q is an N by K orthonormal matrix such that
%  A*Q = Q*S. The diagonal elements of S are eigenvalues of A.
%
%  ... = JDQR('Afun')
%  ... = JDQR('Afun',N)
%  The first input argument is either a square matrix (which can be full
%  or sparse, symmetric or nonsymmetric, real or complex), or a string
%  containing the name of an M-file which applies a linear operator to a
%  given column vector. In the latter case, the M-file must return the 
%  order N of the problem with N = Afun([],'dimension') or N must be
%  specified in the list of input arguments. For example, EIGS('fft',...)
%  is much faster than EIGS(F,...) where F is the explicit FFT matrix.
%
%
%  Other input arguments are optional and can be given in practically any
%  order:
%  ... = JDQR(...,K,SIGMA,OPTIONS,M) 
%  where
%
%      K         An integer, the number of eigenvalues desired.
%      SIGMA     A (vector of) scalar shift(s) or a two letter string.
%      OPTIONS   A structure containing additional parameters.
%      M         A string or a matrix that specifies the preconditioner.
%
%  If K is not specified, then K = MIN(N,5) eigenvalues are computed.
%
%  If SIGMA is not specified, then the K eigenvalues largest in magnitude
%  are computed. If SIGMA is a real or complex scalar, then the K 
%  eigenvalues nearest to SIGMA are computed. If SIGMA is a row of real 
%  or complex scalars, then the I-th eigenvalue that is computed is the 
%  new one that is nearest to SIGMA(1,I1), where I1=MIN(I,LENGTH(SIGMA)).
%  If SIGMA is one of the following strings, then it specifies the desired
%  eigenvalues.
%
%    SIGMA             Location wanted eigenvalues
%
%    'LM'              Largest Magnitude  (the default)
%    'SM'              Smallest Magnitude (same as sigma = 0)
%    'LR'              Largest Real part
%    'SR'              Smallest Real part
%    'BE'              Both Ends.  Computes k/2 eigenvalues
%                      from each end of the spectrum (one more
%                      from the high end if k is odd.)
%
%
%  The OPTIONS structure specifies certain parameters in the algorithm.
%
%    Field name         Parameter                             Default
%
%    OPTIONS.Tol        Convergence tolerance:                 1e-8 
%                          norm(A*Q-Q*S,2) <= tol.   
%    OPTIONS.jmin       Minimum dimension search subspace.     k+5
%    OPTIONS.jmax       Maximum dimension search subspace.     jmin+5
%    OPTIONS.MaxIt      Maximum number of iterations.          100
%    OPTIONS.v0         Starting space.                        ones+0.1*rand
%    OPTIONS.Schur      Gives schur decomposition ('Schur'     'no'
%                         = 'yes') also in case of 2 or 3 
%                         output arguments (X=Q, Lambda=R).
%
%    OPTIONS.TestSpace  For using harmonic Ritz values. If     'Standard'
%                         'TestSpace'='Harmonic' then SIGMA
%                         = 0 is the default value for SIGMA.
%
%    OPTIONS.Disp       If 'Disp'=1 then,                     0
%                           the input parameters, 
%                           the residual size at each step,
%                           and the eigenvalues at detection 
%                         are displayed, and      
%                         the convergence history is plotted.
%                       If 'Disp'=2 then, in addition, 
%                         approximate eigenvalues are plotted 
%                         at each step.
%
%    OPTIONS.LSolver    Linear solver.                         'GMRES'
%    OPTIONS.LS_Tol     Residual reduction linear solver.      1, 0.7, 0.7^2,..
%    OPTIONS.LS_MaxIt   Maximum number it.  linear solver.     5
%    OPTIONS.LS_ell     ell for BiCGstab(ell).                 4
%
%    OPTIONS.Precond    Preconditioner.                        M=[].
%
%  For instance
%
%    OPTIONS=STRUCT('Tol',1.0e-10,'LSolver','BiCGstab',...
%                                           'LS_ell',2,'Precond','M');
%
%  changes the convergence tolerance to 1.0e-10, takes BiCGstab(2) as 
%  linear solver and the preconditioner defined in M.m.
%
%  If M is not specified then there is no preconditioning.
%  The preconditioner can be specified in the OPTIONS structure, but also
%  in the argument list:
%   ... = JDQR(...,K,SIGMA,M,OPTIONS) 
%   ... = JDQR(...,K,SIGMA,L,U,OPTIONS) 
%   ... = JDQR(...,K,SIGMA,'M',OPTIONS)
%   ... = JDQR(...,K,SIGMA,'L','U',OPTIONS)
%  as an N by N matrix M (then M is the preconditioner), or an N by 2*N 
%  matrix M (then  L*U is the preconditioner, where  M = [L,U]), or as
%  N by N matrices L and U (then  L*U is the preconditioner), or as one 
%  or two strings containing the name of  M-files ('M', or 'L' and 'U')
%  which apply a linear operator to a given column vector.
%
%
%  [X,Lambda,HISTORY] = JDQR(A,...) 
%  [X,Lambda,Q,S,HISTORY] = JDQR(A,...) 
%  returns also the convergence history. HISTORY is an array of 3 columns: 
%  HISTORY(I,1) is the residual norm at step J = HISTORY(I,2), HISTORY(I,3)
%  is the cumulative number of multiplications by A at step J. If a search
%  for a new eigenvalue is started at step J then J = HISTORY(I,2) =
%  HISTORY(I+1,2), HISTORY(I,1) is the norm of the "old" residual, and
%  HISTORY(I+1,1) is the norm of the "new" one. HISTORY is empty if the
%  required number of eigenvalues are detected in the initialization phase.
%
%
%  JDQR (without input arguments) lists the options and the defaults.

%  [X,Lambda,Q,S,HISTORY,V] = JDQR(A,...) 
%  returns also the search subspace V at termination.

%   Gerard Sleijpen.
%   Copyright (c) 98

global Qschur Rschur ...
       Operator_MVs Precond_Solves ...
       PinvQ QastPinvQ

if nargin==0, ShowOptions, return, end

%%% Read/set parameters
[n,nselect,sigma,SCHUR,...
   jmin,jmax,tol,maxit,V,INTERIOR,SHOW,PAIRS,JDV0,FIX_tol,t_tol,...
   lsolver,LSpar] = ReadOptions(varargin{1:nargin});
LSpar0=LSpar; JDV=0; tol=tol/sqrt(nselect); tol0=tol; LOCK0=~ischar(sigma); 
if nargout>3 & SCHUR ==1, SCHUR=0; end
tau=0; if INTERIOR>=1 & LOCK0, tau=sigma(1); end
n_tar=size(sigma,1); nt=1; FIG=get(0,'CurrentFigure');
EXPAND=0;

%%% Initiate global variables
Qschur = zeros(n,0); Rschur = []; 
PinvQ  = zeros(n,0); QastPinvQ = [];
Operator_MVs = 0; Precond_Solves = 0; history = [];

%%% Return if eigenvalueproblem is trivial
if n<2
  if n==1, Qschur=1; Rschur=MV(1); end
  if nargout == 0, eigenvalue=Rschur, else
  [varargout{1:nargout}]=output(history,Qschur,Rschur); end, 
return, end

String = ['\r#it=%i #MV=%i dim(V)=%i |r_%i|=%6.1e  '];
StrinP = '--- Checking for conjugate pair ---\n';
time = clock;

%%% Initialize V, W:
%%%   V,W orthonormal, A*V=W*R+Qschur*E, R upper triangular
[V,W,R,E,M]=SetInitialSpaces(V,nselect,tau,jmin,tol); 
j=size(V,2); k=size(Rschur,1);
nit=0; nlit=0;


if k<nselect
switch INTERIOR

case 0

%%% The JD loop (Standard)
%%%    V orthogonal, V orthogonal to Qschur
%%%    V*V=eye(j), Qschur'*V=0, 
%%%    W=A*V, M=V'*W
%%%
W=W*R; if tau ~=0; W=W+tau*V; end, M=M'*R; temptarget=sigma(nt,:);  
while (k<nselect) & (nit < maxit) 

   %%% Compute approximate eigenpair and residual
   [UR,S]=SortSchur(M,temptarget,j==jmax,jmin); 
   y=UR(:,1); theta=S(1,1); u=V*y; w=W*y; 
   r=w-theta*u; [r,s]=RepGS(Qschur,r,0); nr=norm(r); r_KNOWN=1;
   if LOCK0 & nr<t_tol, temptarget=[theta;sigma(nt,:)]; end

          % defekt=abs(norm(RepGS(Qschur,MV(u)-theta*u,0))-nr); 
          % DispResult('defekt',defekt,3)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   history=[history;nr,nit,Operator_MVs];                         %%%
   if SHOW, fprintf(String,nit,Operator_MVs,j,nlit,nr)            %%%
     if SHOW == 2, LOCK =  LOCK0 & nr<t_tol;                       %%%
       if MovieTheta(n,nit,diag(S),jmin,sigma(nt,:),LOCK,j==jmax)  %%%
   break, end, end, end                                            %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   %%% Check for convergence
   if nr<tol

      %%% Expand the partial Schur form
      Qschur=[Qschur,u]; 
      %% Rschur=[[Rschur;zeros(1,k)],Qschur'*MV(u)]; k=k+1; 
      Rschur=[Rschur,s;zeros(1,k),theta];  k=k+1;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if SHOW, ShowLambda(theta,k), end %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if k>=nselect, break, end, r_KNOWN=0; JDV=0;
      

      %%% Expand preconditioned Schur matrix PinvQ
      UpdateMinv;

      if j==1, 
         [V,W,R,E,M]=SetInitialSpaces(zeros(n,0),nselect,tau,jmin,tol); 
         k=size(Rschur,1); if k>=nselect, break, end
         W=W*R; if tau ~=0; W=W+tau*V; end; M=M'*R; j=size(V,2);
      else,
         J=[2:j]; j=j-1; UR=UR(:,J); 
         M=S(J,J); V=V*UR; W=W*UR; 
      end

     if PAIRS & abs(imag(theta))>tol, v=imag(u/sign(max(u)));
       if norm(v)>tol, v=RepGS(Qschur,v,0); EXPAND=(norm(v)>sqrt(tol)); end
     end

     if EXPAND, temptarget=conj(theta); if SHOW, fprintf(StrinP), end
     else, nlit=0; nt=min(nt+1,n_tar); temptarget=sigma(nt,:); end

   end % nr<tol

   %%% Check for shrinking the search subspace
   if j>=jmax
      j=jmin; J=[1:j]; UR=UR(:,J);
      M=S(J,J); V=V*UR; W=W*UR;
   end % if j>=jmax
 
   if r_KNOWN
      if JDV, disp('Stagnation'), u=V; end
      %%% Solve correction equation
      if FIX_tol*nr>1 & LOCK0, theta=tau; else, FIX_tol=0; end
      v=Solve_pce(theta,u,r,lsolver,LSpar,nlit);
      nlit=nlit+1; nit=nit+1; r_KNOWN=0; EXPAND=1; JDV=0;
   end % if r_KNOWN

   if EXPAND
     %%% Expand the subspaces of the interaction matrix  
     [v,zeta]=RepGS([Qschur,V],v);
     if JDV0 & abs(zeta(end,1))/norm(zeta)<0.06, JDV=JDV+1; end
     if size(v,2)>0
        w=MV(v);
        M=[M,V'*w;v'*W,v'*w]; 
        V=[V,v]; W=[W,w]; j=j+1; EXPAND=0; tol=tol0;
      else
        tol=2*tol;
      end
   end % if EXPAND

end % while (nit<maxit)


case 1

%%% The JD loop (Harmonic Ritz values)
%%%    Both V and W orthonormal and orthogonal w.r.t. Qschur
%%%    V*V=eye(j), Qschur'*V=0, W'*W=eye(j), Qschur'*W=0
%%%    (A*V-tau*V)=W*R+Qschur*E, E=Qschur'*(A*V-tau*V), M=W'*V
%%%
temptarget=0; lsolver0=lsolver;
while (k<nselect) & (nit<maxit) 

   %%% Compute approximate eigenpair and residual
   [UR,UL,S,T]=SortQZ(R,M,temptarget,j>=jmax,jmin);
   y=UR(:,1); theta=T(1,1)'*S(1,1); 
   u=V*y; w=W*(R*y); r=w-theta*u; nr=norm(r); r_KNOWN=1; 
   if nr<t_tol, temptarget=[theta;0]; end, theta=theta+tau;

           % defekt=abs(norm(RepGS(Qschur,MV(u)-theta*u,0))-nr); 
           % DispResult('defect',defekt,3)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   history=[history;nr,nit,Operator_MVs];                           %%%
   if SHOW, fprintf(String,nit,Operator_MVs,j,nlit,nr)              %%%
     if SHOW == 2, Lambda=diag(S)./diag(T)+tau; Lambda(1)=theta;     %%%
       if MovieTheta(n,nit,Lambda,jmin,sigma(nt,:),nr<t_tol,j==jmax) %%%
   break, end, end, end                                              %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%% Check for convergence
   if nr<tol

     %%% Expand the partial Schur form
     Qschur=[Qschur,u]; 
     %% Rschur=[[Rschur;zeros(1,k)],Qschur'*MV(u)]; k=k+1;
     Rschur=[Rschur,E*y;zeros(1,k),theta];   k=k+1;  

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if SHOW, ShowLambda(theta,k), end %%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if k>=nselect, break, end, r_KNOWN=0; JDV=0;

     %%% Expand preconditioned Schur matrix PinvQ
     UpdateMinv;

     if j==1,
       [V,W,R,E,M]=...
         SetInitialSpaces(zeros(n,0),nselect,tau,jmin,tol); 
       k=size(Rschur,1); if k>=nselect, break, end, j=size(V,2);
     else
       J=[2:j]; j=j-1; UR=UR(:,J); UL=UL(:,J);
       R=S(J,J); M=T(J,J); V=V*UR; W=W*UL;    
       [r,a]=RepGS(u,r,0);      E=[E*UR;(T(1,1)'-a/S(1,1))*S(1,J)];
                                     
       s=(S(1,J)/S(1,1))/R; W=W+r*s; M=M+s'*(r'*V);
       if (nr*norm(s))^2>eps, [W,R0]=qr(W,0); R=R0*R; M=R0'\M; end

     end

     if PAIRS & abs(imag(theta))>tol, v=imag(u/sign(max(u)));
       if norm(v)>tol, v=RepGS(Qschur,v,0); EXPAND=(norm(v)>sqrt(tol)); end
     end

     if EXPAND, if SHOW, fprintf(StrinP), end
       temptarget=[conj(theta)-tau;0];
     else, nlit=0; temptarget=0;
       if nt<n_tar
         nt=nt+1; tau0=tau; tau=sigma(nt,1); tau0=tau0-tau;
         [W,R]=qr(W*R+tau0*V,0); M=W'*V;
       end
     end

   end

   %%% Check for shrinking the search subspace
   if j>=jmax
      j=jmin; J=[1:j]; UR=UR(:,J); UL=UL(:,J);
      R=S(J,J); M=T(J,J); V=V*UR; W=W*UL;           E=E*UR;
   end % if j>=jmax

   if r_KNOWN
      %%% Solve correction equation
      if JDV, disp('Stagnation'), u=V;
        % LSpar(end-1)=(LSpar(end-1)+15)*2;
        % lsolver='bicgstab'; LSpar=[1.e-2,300,4];
      else
        % LSpar=LSpar0; JDV=0; lsolver=lsolver0;
      end
      if FIX_tol*nr>1 & LOCK0, theta=tau; else, FIX_tol=0; end
      v=Solve_pce(theta,u,r,lsolver,LSpar,nlit);
      nlit=nlit+1; nit=nit+1; r_KNOWN=0; EXPAND=1; JDV=0;
   end
 
   if EXPAND
      %%% Expand the subspaces of the interaction matrix  
      [v,zeta]=RepGS([Qschur,V],v);
      if JDV0 & abs(zeta(end,1))/norm(zeta)<0.06, JDV=JDV+1; end
      if size(v,2)>0 
        w=MV(v); if tau ~=0, w=w-tau*v; end
        [w,e]=RepGS(Qschur,w,0); [w,y]=RepGS(W,w); 
        R=[[R;zeros(1,j)],y]; M=[M,W'*v;w'*V,w'*v];  E=[E,e];
        V=[V,v]; W=[W,w]; j=j+1; EXPAND=0; tol=tol0;
      else
        tol=2*tol;
      end
   end
     
end % while (nit<maxit)

case 1.1

%%% The JD loop (Harmonic Ritz values)
%%%    V W AV.
%%%    Both V and W orthonormal and orthogonal w.r.t. Qschur, AV=A*V-tau*V
%%%    V*V=eye(j),  W'*W=eye(j), Qschur'*V=0, Qschur'*W=0, 
%%%    (I-Qschur*Qschur')*AV=W*R, M=W'*V; R=W'*AV;
%%%
AV=W*R; temptarget=0;
while (k<nselect) & (nit<maxit) 

   %%% Compute approximate eigenpair and residual
   [UR,UL,S,T]=SortQZ(R,M,temptarget,j>=jmax,jmin);
   y=UR(:,1); u=V*y; w=AV*y; theta=u'*w;  
   r=w-theta*u; [r,y]=RepGS(Qschur,r,0); nr=norm(r); r_KNOWN=1;
   if nr<t_tol, temptarget=[theta;0]; end, theta=theta+tau;

           % defekt=abs(norm(RepGS(Qschur,MV(u)-theta*u,0))-nr); 
           % DispResult('defect',defekt,3)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   history=[history;nr,nit,Operator_MVs];                           %%%
   if SHOW, fprintf(String,nit,Operator_MVs,j,nlit,nr)              %%%
     if SHOW == 2,Lambda=diag(S)./diag(T)+tau; Lambda(1)=theta;      %%%
       if MovieTheta(n,nit,Lambda,jmin,sigma(nt,:),nr<t_tol,j==jmax) %%%
   break, end, end, end                                              %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   %%% Check for convergence
   if nr<tol

     %%% Expand the partial Schur form
     Qschur=[Qschur,u]; 
     %% Rschur=[[Rschur;zeros(1,k)],Qschur'*MV(u)]; k=k+1;
     Rschur=[Rschur,y;zeros(1,k),theta];  k=k+1; 

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if SHOW, ShowLambda(theta,k), end %%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if k>=nselect, break, end, r_KNOWN=0;

     %%% Expand preconditioned Schur matrix PinvQ
     UpdateMinv;

     if j==1
       [V,W,R,E,M]=SetInitialSpaces(zeros(n,0),nselect,tau,jmin,tol); 
       k=size(Rschur,1); if k>=nselect, break, end
       AV=W*R; j=size(V,2); 
     else
       J=[2:j]; j=j-1; UR=UR(:,J); UL=UL(:,J);
       AV=AV*UR; R=S(J,J); M=T(J,J); V=V*UR; W=W*UL;
     end
      
     if PAIRS & abs(imag(theta))>tol, v=imag(u/sign(max(u)));
       if norm(v)>tol, v=RepGS(Qschur,v,0); EXPAND=(norm(v)>sqrt(tol)); end
     end

     if EXPAND,if SHOW, fprintf(StrinP), end
       temptarget=[conj(theta)-tau;0];
     else, nlit=0; temptarget=0;
       if nt<n_tar
         nt=nt+1; tau0=tau; tau=sigma(nt,1); tau0=tau0-tau;
         AV=AV+tau0*V; [W,R]=qr(W*R+tau0*V,0); M=W'*V;
       end
     end

   end

   %%% Check for shrinking the search subspace
   if j>=jmax
     j=jmin; J=[1:j]; UR=UR(:,J); UL=UL(:,J);
     AV=AV*UR; R=S(J,J); M=T(J,J); V=V*UR; W=W*UL;
   end % if j>=jmax

   if r_KNOWN
     %%% Solve correction equation
     if FIX_tol*nr>1 & LOCK0, theta=tau; else, FIX_tol=0; end
     v=Solve_pce(theta,u,r,lsolver,LSpar,nlit);
     nlit=nlit+1; nit=nit+1; r_KNOWN=0; EXPAND=1;
   end

   if EXPAND
     %%% Expand the subspaces of the interaction matrix  
     v=RepGS([Qschur,V],v); 
     if size(v,2)>0
       w=MV(v); if tau ~=0, w=w-tau*v;end
       AV=[AV,w]; R=[R,W'*w];
       w=RepGS([Qschur,W],w);
       R=[R;w'*AV]; M=[M,W'*v;w'*V,w'*v]; 
       V=[V,v]; W=[W,w]; j=j+1; EXPAND=0; tol=tol0;
     else
       tol=2*tol;
     end
   end
     
end % while (nit<maxit)

case 1.2

%%% The JD loop (Harmonic Ritz values)
%%%    W orthonormal, V and W orthogonal to Qschur, 
%%%    W'*W=eye(j), Qschur'*V=0, Qschur'*W=0
%%%    W=(A*V-tau*V)-Qschur*E, E=Qschur'*(A*V-tau*V), 
%%%    M=W'*V
V=V/R; M=M/R; temptarget='LM';            E=E/R;
while (k<nselect) & (nit<maxit) 

   %%% Compute approximate eigenpair and residual
   [UR,S]=SortSchur(M,temptarget,j==jmax,jmin);
   y=UR(:,1); u=V*y; nrm=norm(u); y=y/nrm; u=u/nrm;
   theta=S(1,1)'/(nrm*nrm); w=W*y; r=w-theta*u; nr=norm(r); r_KNOWN=1;
   if nr<t_tol, temptarget=[S(1,1);inf]; end, theta=theta+tau;

           % defekt=abs(norm(RepGS(Qschur,MV(u)-theta*u,0))-nr); 
           % DispResult('defect',defekt,3)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   history=[history;nr,nit,Operator_MVs];                           %%%
   if SHOW, fprintf(String,nit,Operator_MVs,j,nlit,nr)              %%%
     if SHOW == 2, Lambda=1./diag(S)+tau; Lambda(1)=theta;           %%%
       if MovieTheta(n,nit,Lambda,jmin,sigma(nt,:),nr<t_tol,j==jmax) %%%
   break, end, end, end                                              %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%% Check for convergence
   if nr<tol

     %%% Expand the partial Schur form
     Qschur=[Qschur,u]; 
     %% Rschur=[[Rschur;zeros(1,k)],Qschur'*MV(u)]; k=k+1;
     y=E*y; Rschur=[Rschur,y;zeros(1,k),theta]; k=k+1;  

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if SHOW, ShowLambda(theta,k), end %%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if k>=nselect, break, end, r_KNOWN=0;

     %%% Expand preconditioned Schur matrix PinvQ
     UpdateMinv;

     if j==1
       [V,W,R,E,M]=SetInitialSpaces(zeros(n,0),nselect,tau,jmin,tol); 
       k=size(Rschur,1); if k>=nselect, break, end
       V=V/R; j=size(V,2);  M=M/R;         E=E/R;
     else
       J=[2:j]; j=j-1; UR=UR(:,J); M=S(J,J);
       V=V*UR; W=W*UR; [r,a]=RepGS(u,r,0); 
       s=u'*V; V=V-u*s; W=W-r*s; M=M-s'*(r'*V)-(W'*u)*s;     
                                 E=[E*UR-y*s;(tau-theta-a)*s];
         
       if (nr*norm(s))^2>eps, [W,R]=qr(W,0); V=V/R; M=(R'\M)/R; E=E/R; end

     end

     if PAIRS & abs(imag(theta))>tol, v=imag(u/sign(max(u)));
       if norm(v)>tol, v=RepGS(Qschur,v,0); EXPAND=(norm(v)>sqrt(tol)); end
     end

     if EXPAND, if SHOW, fprintf(StrinP), end
       temptarget=[1/(conj(theta)-tau);inf];
     else, nlit=0; temptarget='LM';
       if nt<n_tar
         nt=nt+1; tau0=tau; tau=sigma(nt,1); 
         [W,R]=qr(W+(tau0-tau)*V,0); V=V/R; M=W'*V; E=E/R; 
       end
     end

   end

   %%% Check for shrinking the search subspace
   if j>=jmax
     j=jmin; J=[1:j]; UR=UR(:,J);
     M=S(J,J); V=V*UR; W=W*UR;               E=E*UR;
   end % if j>=jmax

   if r_KNOWN
     %%% Solve correction equation
     if FIX_tol*nr>1 & LOCK0, theta=tau; else, FIX_tol=0; end
     v=Solve_pce(theta,u,r,lsolver,LSpar,nlit);
     nlit=nlit+1; nit=nit+1; r_KNOWN=0; EXPAND=1;
   end

   if EXPAND
     %%% Expand the subspaces of the interaction matrix  
     v=RepGS(Qschur,v,0);
     if size(v,2)>0 
       w=MV(v); if tau ~=0, w=w-tau*v; end
       [w,e]=RepGS(Qschur,w,0); [w,y]=RepGS(W,w); 
       nrw=y(j+1,1); y=y(1:j,:);   
       v=v-V*y; v=v/nrw;                     e=e-E*y; e=e/nrw;
       M=[M,W'*v;w'*V,w'*v]; 
       V=[V,v]; W=[W,w]; j=j+1;              E=[E,e];

       if 1/cond(M)<10*tol
         [V,W,R,E,M]=SetInitialSpaces(V,nselect,tau,jmin,tol,W,E);    
         k=size(Rschur,1); if k>=nselect, break, end
         V=V/R; M=M/R; j=size(V,2);temptarget='LM'; E=E/R;
       end 

       EXPAND=0; tol=tol0;
     else
       tol=2*tol;
     end
   end 

end % while (nit<maxit)

end % case
end % if (k<nselect)

time_needed=etime(clock,time);

Refine([Qschur,V],1);% 2-SCHUR);
CheckSortSchur(sigma);

Lambda=[]; X=zeros(n,0); 
if SCHUR ~= 1 & k>0, [z,Lambda]=Jordan(Rschur,SCHUR); X=Qschur*z; end

%-------------- display results ----------------------------
if SHOW == 2 & nit>0, 
   MovieTheta, if isempty(FIG), FIG=1; end, figure(FIG), 
end
if SHOW & size(history,1)>0

   switch INTERIOR
      case 0
        testspace='V, V orthonormal';     
      case 1
        testspace='A*V-sigma*V, V and W orthonormal';    
      case 1.1
        testspace='A*V-sigma*V, V and W orthonormal, AV';   
      case 1.2
        testspace='A*V-sigma*V, W orthogonal';     
      otherwise 
        testspace='Experimental'; 
   end

   StringT=sprintf('The test subspace W is computed as  W = %s.',testspace);
   StringX=sprintf('JDQZ with jmin=%g, jmax=%g, residual tolerance %g.',...
           jmin,jmax,tol); 
   StringY=sprintf('Correction equation solved with %s.',lsolver);   
   date=fix(clock);
   String=sprintf('\n%2i-%2i-%2i, %2i:%2i:%2i',date(3:-1:1),date(4:6));

   StringL='log_{10} || r_{#it} ||_2'; if SHOW, SHOW=2; end
   for pl=1:SHOW
     subplot(SHOW,1,pl), t=history(:,pl+1);
     plot(t,log10(history(:,1)),'*-',t,log10(tol)+0*t,':')
     legend(StringL), title(StringT)
     StringL='log_{10} || r_{#MV} ||_2'; StringT=StringX; 
   end 
   if SHOW==2, xlabel([StringY,String])
   else,  xlabel([StringX,String]), ylabel(StringY), end
   drawnow
end

if SHOW
   str1=num2str(abs(k-nselect)); str='s';
   if k>nselect, 
     if k==nselect+1, str1='one'; str=''; end
     fprintf('\n\nDetected %s additional eigenpair%s.',str1,str)
   end
   if k<nselect, 
     if k==0, str1='any'; str=''; elseif k==nselect-1, str1='one'; str=''; end
     fprintf('\n\nFailed detection of %s eigenpair%s.',str1,str)
   end
   if k>0, ShowLambda(diag(Rschur)); else, fprintf('\n'); end

   Str='time_needed';                              DispResult(Str,eval(Str))
   fprintf('\n%36s: %i','Number of Operator actions',Operator_MVs)
   if Precond_Solves
   fprintf('\n%36s: %i','Number of preconditioner solves',ceil(Precond_Solves))
   end

   if (k>0)
      if SCHUR ~= 1
        Str='norm(MV(X)-X*Lambda)';                DispResult(Str,eval(Str))
      end
      Str='norm(MV(Qschur)-Qschur*Rschur)';        DispResult(Str,eval(Str))
      I=eye(k); Str='norm(Qschur''*Qschur-I)';     DispResult(Str,eval(Str))  
   end
   fprintf('\n\n')

end

na=nargout; if na == 6, varargout{6}=V; na=5; end 
if ~na, if ~SHOW, eigenvalues=diag(Rschur), end, return, end
[varargout{1:na}]=output(history,X,Lambda);

return

%===========================================================================
%======= PREPROCESSING =====================================================
%===========================================================================

%======= INITIALIZE SUBSPACE ===============================================
function [V,W,R,E,M]=SetInitialSpaces(V,nselect,tau,jmin,tol,W,E);
%[V,W,R,E,M]=SetInitialSpaces(VV,nselect,tau,jmin,tol);
%  Output: V(:,1:SIZE(VV,2))=ORTH(VV),
%          V'*V=W'*W=EYE(JMIN), M=W'*V;
%          such that A*V-tau*V=W*R+Qschur*E, 
%          with R upper triangular, and E=Qschur'*(A*V-tau*V).
%
%[V,W,R,E,M]=SetInitialSpaces(VV,nselect,tau,jmin,tol,AV,EE);
%  Input such that
%  A*VV-tau*VV=AV+Qschur*EE, EE=Qschur'*(A*VV-tau*VV);
%
%  Output: V(:,1:SIZE(VV,2))=ORTH(VV),
%          V'*V=W'*W=EYE(JMIN), M=W'*V;
%          such that A*V-tau*V=W*R+Qschur*E, 
%          with R upper triangular, and E=Qschur'*(A*V-tau*V).

global Qschur Rschur

[n,j]=size(V); k=size(Qschur,2);

if j>1, 
  [V,R]=qr(V,0);
  if nargin <6,   
    W=MV(V); R=eye(j); if tau~=0, W=W-tau*V; end
    if k>0, E=Qschur'*W; W=W-Qschur*E; else, E=zeros(0,j); end
  end
  [V,W,R,E,M]=CheckForNullSpace(V,nselect,tau,tol,W,E,R);
  l=size(Qschur,2); j=size(V,2);
  if l>=nselect, if size(V,2)==0; R=1; M=1; return, end, end
  if l>k, UpdateMinv; end, k=l;
end
if j==0, nr=0;
  while nr==0
    V = ones(n,1)+0.1*rand(n,1); V=RepGS(Qschur,V); nr=norm(V);
  end, j=1;
end
if j==1
  [V,H,E]=Arnoldi(V,tau,jmin,nselect,tol); 
  l=size(Qschur,2); j=max(size(H,2),1);
  if l>=nselect, W=V; R=eye(j+1); M=R; return, end
  if l>k, UpdateMinv; end
  [Q,R]=qr(full(H),0);
  W=V*Q; V(:,j+1)=[]; M=Q(1:j,:)';
  %% W=V*Q; V=V(:,1:j)/R; E=E/R; R=eye(j); M=Q(1:j,:)'/R;
  %% W=V*H; V(:,j+1)=[];R=R'*R;   M=H(1:j,:)';
end

return
%%%======== ARNOLDI (for initializing spaces) ===============================
function [V,H,E]=Arnoldi(v,tau,jmin,nselect,tol)
%
%[V,AV,H,nMV,tau]=ARNOLDI(A,V0,TAU,JMIN,NSELECT,TOL)
%    ARNOLDI computes the Arnoldi factorization of dimenison JMIN+1:
%    (A-tau)*V(:,1:JMIN)=V*H where V is n by JMIN+1 orthonormal with
%    first column a multiple of V0, and H is JMIN+1 by JMIN Hessenberg.
%
%    If an eigenvalue if H(1:j,1:j) is an eigenvalue of A
%    within the required tolerance TOL then the Schurform
%    A*Qschur=Qschur*Rschur is expanded and the Arnoldi factorization
%    (A-tau)*V(:,1:j)=V(:,1:j+1)*H(1:j+1,1:j) is deflated.
%    Returns if size(Qschur,2) = NSELECT or size(V,2) = JMIN+1

%    (A-tau)*V(:,1:JMIN)=V*H+Qschur*E, Qschur'*V=0

%   Coded November 5, 1998, G. Sleijpen

global Qschur Rschur

k=size(Qschur,2); [n,j]=size(v);

if ischar(tau), tau=0; end

H=zeros(1,0); V=zeros(n,0); E=[];
j=0; nr=norm(v);

while j<jmin & k<nselect & j+k<n
   if nr>=tol
      v=v/nr; V=[V,v]; j=j+1;
      Av=MV(v);
   end
   if j==0 
      H=zeros(1,0); j=1;
      nr=0; while nr==0, v=RepGS(Qschur,rand(n,1)); nr=norm(v); end
      v=v/nr; V=v; Av=MV(v);
   end
   if tau~=0; Av=Av-tau*v; end, [v,e] = RepGS(Qschur,Av,0);
   if k==0, E=zeros(0,j); else, E = [E,e(1:k,1)]; end
   [v,y] = RepGS(V,v,0);        H = [H,y(1:j,1)];
   nr = norm(v);                H = [H;zeros(1,j-1),nr];
   [Q,U,H1] = DeflateHess(full(H),tol); 
   j=size(U,2); l=size(Q,2);
   if l>0 %--- expand Schur form ------
      Qschur=[Qschur,V*Q]; 
      Rschur=[Rschur,E*Q; zeros(l,k),H1(1:l,1:l)+tau*eye(l)]; k=k+l;
      E=[E*U;H1(1:l,l+1:l+j)];
      if j>0, V=V*U; H=H1(l+1:l+j+1,l+1:l+j);
      else, V=zeros(n,0); H=zeros(1,0); end
   end
end % while

if nr>=tol
   v=v/nr; V=[V,v]; 
end

return
%----------------------------------------------------------------------
function [Q,U,H]=DeflateHess(H,tol)
% H_in*[Q,U]=[Q,U]*H_out such that H_out(K,K) upper triangular
% where K=1:SIZE(Q,2) and ABS(Q(end,2)*H_in(j+1,j))<TOL, j=SIZE(H,2),

[j1,j]=size(H); 
if j1==j, [Q,H]=schur(H); U=zeros(j,0); return, end

nr=H(j+1,j);
U=eye(j); i=1; J=i:j;
for l=1:j
   [X,Lambda]=eig(H(J,J));
   I=find(abs(X(size(X,1),:)*nr)<tol);
   if isempty(I), break, end
   q=X(:,I(1)); q=q/norm(q); 
   q(1,1)=q(1,1)+sign(q(1,1)); q=q/norm(q);
   H(:,J)=H(:,J)-(2*H(:,J)*q)*q'; 
   H(J,:)=H(J,:)-2*q*(q'*H(J,:));
   U(:,J)=U(:,J)-(2*U(:,J)*q)*q'; 
   i=i+1; J=i:j;
end

[Q,HH]=RestoreHess(H(i:j+1,J)); 
H(:,J)=H(:,J)*Q; H(J,:)=Q'*H(J,:);
U(:,J)=U(:,J)*Q; 

Q=U(:,1:i-1); U=U(:,i:j);

return
%----------------------------------------------------------------------
function [Q,M]=RestoreHess(M)

[j1,j2]=size(M); Q=eye(j2);

for j=j1:-1:2
   J=1:j-1;
   q=M(j,J)'; q=q/norm(q);
   q(j-1,1)=q(j-1,1)+sign(q(j-1,1));
   q=q/norm(q);
   M(:,J)=M(:,J)-2*(M(:,J)*q)*q';
   M(J,:)=M(J,:)-2*q*q'*M(J,:);
   Q(:,J)=Q(:,J)-2*Q(:,J)*q*q';
end

return
%%%=========== END ARNOLDI ============================================  
function [V,W,R,E,M]=CheckForNullSpace(V,nselect,tau,tol,W,E,Rv);
% V,W orthonormal, A*V-tau*V=W*R+Qschur'*E

global Qschur Rschur

  k=size(Rschur,1); j=size(V,2); 
 
  [W,R]=qr(W,0); E=E/Rv; R=R/Rv; M=W'*V; 
  %%% not accurate enough M=Rw'\(M/Rv);

  if k>=nselect, return, end

  CHECK=1; l=k;
  mp=10*sqrt(size(V,1))*eps;

  [S,T,Z,Q]=qz(R,M); Z=Z';
  while CHECK 
    I=SortEigPairVar(S,T,2); [Q,Z,S,T]=SwapQZ(Q,Z,S,T,I(1)); 
    s=abs(S(1,1)); t=min(abs(T(1,1)),1);  % CHECK=(s*sqrt(1-t*t)<tol);
    CHECK=(s*sqrt(1-t*t)<max(tol,sqrt(mp*t*s)));
    % resn1=[s*sqrt(1-t*t),tol,sqrt(mp*t*s)]
    if CHECK 
      V=V*Q; W=W*Z; E=E*Q; 
      u=V(:,1); [r,a]=RepGS(u,(W(:,1)-T(1,1)'*u)*S(1,1),0);
      CHECK=(norm(r)<tol);
      if CHECK
        Qschur=[Qschur,u]; t=(T(1,1)'-a/S(1,1))*S(1,:);
        Rschur=[Rschur,E(:,1);zeros(1,k),tau+t(1,1)]; k=k+1;
        J=[2:j]; j=j-1; 
        V=V(:,J); W=W(:,J); E=[E(:,J);t(1,J)];  s=S(1,J)/S(1,1);
        R=S(J,J); M=T(J,J); Q=eye(j); Z=eye(j); 

        s=s/R; nrs=norm(r)*norm(s);
        if nrs>=tol,
          W=W+r*s; M=M+s'*(r'*V); 
          if nrs^2>eps, [W,R0]=qr(W,0); R=R0*R; M=R0'\M; end
        end

        S=R; T=M;
        CHECK=(k<nselect & j>0);
      end
    end
  end

return


%===========================================================================
%======= POSTPROCESSING ====================================================
%===========================================================================
function Refine(V,gamma);

if gamma==0, return, end

global Qschur Rschur

  J=1:size(Rschur,1); 

  if gamma==1, 
    [V,R]=qr(V(:,J),0); W=MV(V); M=V'*W;
    [U,Rschur]=schur(M);
    [U,Rschur]=rsf2csf(U,Rschur); Qschur=V*U;
    return 
  elseif gamma==2
    [V,R]=qr(V,0); W=MV(V); M=V'*W;
    [U,S]=schur(M); [U,S]=rsf2csf(U,S); 
    R=R*U; F=R'*R-S'*S;
    [X,Lambda]=Jordan(S,0); 
    % Xinv=inv(X); D=sqrt(diag(Xinv*Xinv')); X=X*diag(D);
    [d,I]=sort(abs(diag(X'*F*X)));
    [U,S]=SwapSchur(U,S,I(J));
    Qschur=V*U(:,J); Rschur=S(J,J);
  end

return

%===========================================================================
function CheckSortSchur(sigma)
global Qschur Rschur

k=size(Rschur,1); if k==0, return, end

I=SortEig(diag(Rschur),sigma);

if ~min((1:k)'==I)
   [U,Rschur]=SwapSchur(eye(k),Rschur,I);
   Qschur=Qschur*U;
end

return

%%%=========== COMPUTE SORTED JORDAN FORM ==================================
function [X,Jor]=Jordan(S,SCHUR)
% [X,J]=FINDJORDAN(S)
%   For S k by k upper triangular matrices with ordered diagonal
%   elements.
%   FINDJORDAN computes the Jordan decomposition S*X=X*J, where
%   X is a k by k matrix of eigenvectors and principal vectors,
%   J is a k by k matrix, J is Jordan, i.e. only non-zeros
%   on the diagonal and the first upper diagonal.
%   diag(J) are the eigenvalues.
%   The columns of X are normalized.

% coded by Gerard Sleijpen, Januari 14, 1998

k=size(S,1); 
if k<1, 
  if k==0, X=[];Jor=[]; end
  if k==1, X=1; Jor=S; end
  return
end

tol=k*norm(S,1)*eps;
[X,Jor,I]=PseudoJordan(S,tol);

if SCHUR == 0
for l=1:length(I)-1
  if I(l)<I(l+1)-1, 
    J=[I(l):I(l+1)-1];  
    [U,JJor]=JordanBlock(Jor(J,J),tol);
    X(:,J)=X(:,J)*U; Jor(J,J)=JJor;
  end
end
end

Jor=Jor+diag(diag(S)); 
Jor=Jor.*(abs(Jor)>tol);

return
%==================================================
function [X,Jor,J]=PseudoJordan(S,delta);
% Computes a pseudo-Jordan decomposition for
% the upper triangular matrix S with ordered diagonal elements.
% S*X=X*(diag(diag(S))+Jor) with X(:,i:j) orthonormal if  
% its columns span an invariant subspace of S.

k=size(S,1); s=diag(S); T=eye(k);
Jor=zeros(k); X=eye(k); J=1;
for i=2:k
  I=[1:i]; 
  C=S(I,I)-s(i,1)*T(I,I); 
  C(i,i)=norm(C,inf);
  if C(i,i)>0
    tol=delta*C(i,i);
    for j=i:-1:1 
      if j==1 | abs(C(j-1,j-1))>tol, break; end
    end
    e=zeros(i,1); e(i,1)=1;
    if j==i
      J=[J,i]; q=C\e; X(I,i)=q/norm(q);
    else
      q=X(I,j:i-1); 
      q=[C,T(I,I)*q;q',zeros(i-j)]\[e;zeros(i-j,1)];
      q=q/norm(q(I,1)); X(I,i)=q(I,1);
      Jor(j:i-1,i)=-q(i+1:2*i-j,1);
    end
  end
end
J=[J,k+1];
return
%==================================================
function [X,Jor,U]=JordanBlock(A,tol)
%  If A is nilpotent, then A*X=X*Jor with
%  Jor a Jordan block
%

k=size(A,1); Id=eye(k); 
U=Id; aa=A; j=k; jj=[]; J=1:k;

while j>0
  [u,s,v]=svd(aa); U(:,J)=U(:,J)*v;
  sigma=diag(s); delta=tol;
  J=find(sigma<delta); 
  if isempty(J),j=0; else, j=min(J)-1; end
  jj=[jj,j]; if j==0, break, end
  aa=v'*u*s; J=1:j; aa=aa(J,J);
end 
Jor=U'*A*U;
Jor=Jor.*(abs(Jor)>tol);
l=length(jj); jj=[jj(l:-1:1),k];

l2=jj(2)-jj(1); J=jj(1)+(1:l2); 
JX=Id(:,J); X=Id;
for j=2:l
  l1=l2+1; l2=jj(j+1)-jj(j); 
  J2=l1:l2; J=jj(j)+(1:l2);
  JX=Jor*JX; D=diag(sqrt(diag(JX'*JX))); JX=JX/D;
  [Q,S,V]=svd(JX(J,:));
  JX=[JX,Id(:,J)*Q(:,J2)]; X(:,J)=JX;
end

J=[];
for i=1:l2
  for k=l:-1:1
    j=jj(k)+i; if j<=jj(k+1), J=[J,j]; end
  end
end

X=X(:,J); Jor=X\(Jor*X); X=U*X;
Jor=Jor.*(abs(Jor)>100*tol);

return
%========== OUTPUT =========================================================
function varargout=output(history,X,Lambda)

global Qschur Rschur

if nargout == 1, varargout{1}=diag(Rschur);         return, end
if nargout == 3, varargout{3}=history;                      end
if nargout > 3,  varargout{3}=Qschur;  varargout{4}=Rschur; end
if nargout > 4,  varargout{5}=history;                      end
if nargout < 4 & size(X,2)<2
  varargout{1}=Qschur; varargout{2}=Rschur; return
else
  varargout{1}=X; varargout{2}=Lambda;
end

return

%===========================================================================
%===== UPDATE PRECONDITIONED SCHUR VECTORS =================================
%===========================================================================
function   UpdateMinv

global Qschur PinvQ QastPinvQ Precond_L

   if ~isempty(Precond_L)
     [n,k]=size(Qschur); l=size(PinvQ,2);
     if l==0, PinvQ=zeros(n,0); QastMinvZ = []; end
     Pinv_u=SolvePrecond(Qschur(:,l+1:k)); 
     QastPinvQ=[[QastPinvQ;Qschur(:,k)'*PinvQ],Qschur'*Pinv_u];
     PinvQ=[PinvQ,Pinv_u]; 
   end

return
%===========================================================================
%===== SOLVE CORRECTION EQUATION ===========================================
%===========================================================================
function t=Solve_pce(theta,u,r,lsolver,par,nit)

global Qschur PinvQ QastPinvQ Precond_L

  Q=[Qschur,u];

  switch lsolver
    case {'exact','iluexact'}

      t = feval(lsolver,theta,Q,r);

    case {'gmres','bicgstab','olsen','cg','minres','symmlq'}
  
         %%% compute vectors and matrices for skew projection
         if isempty(Precond_L)
           PQ=Q; QPQ=1;
         else
           [PQ,QPQ]=FormPQ(u);
         end

         %%% solve preconditioned system
         t = feval(lsolver,theta,Q,PQ,QPQ,r,spar(par,nit));
      
  end
    
return
%------------------------------------------------------------------------
function [PQ,QPQ]=FormPQ(q)
% compute vectors and matrices for skew projection

global Qschur PinvQ QastPinvQ Precond_L

if isempty(Precond_L)
  PQ=[Qschur,q]; QPQ=1;
else
  Pinv_q=SolvePrecond(q);
  QPQ=[QastPinvQ,Qschur'*Pinv_q;q'*PinvQ,q'*Pinv_q];
  PQ=[PinvQ,Pinv_q];
end
    
return
%=======================================================================
%======= LINEAR SOLVERS ================================================
%=======================================================================
function x = exact(theta,Q,r)

global Operator_A

   [n,k]=size(Q); [n,l]=size(r);

   if ischar(Operator_A)
      [x,xtol]=bicgstab(theta,Q,Q,1,r,[5.0e-14/norm(r),200,4]);
      return
   end
   
   Aug=[Operator_A-theta*speye(n,n),Q;Q',zeros(k,k)];
   x = Aug\[r;zeros(k,l)]; x = x(1:n,1:l);  

return
%----------------------------------------------------------------------
function x = iluexact(theta,Q,r)

global Precond_L Precond_U

   [n,k]=size(Q); [n,l]=size(r);

   y = Precond_L\[r,Q]; 
   x = [Precond_U,y(:,l+1:k+1);Q',zeros(k,k)]\[y(:,1:l);zeros(k,l)]; 
   x = x(1:n,1:l); 
 
return
%----------------------------------------------------------------------
function r = olsen(theta,Q,PQ,QPQ,r,par)
% returns the preconditioned residual as approximate solution
% May be sufficient in case of an excellent preconditioner
   r=SkewProj(Q,PQ,QPQ,SolvePrecond(r));
return
%
%======= Iterative methods =============================================
%------------------------------------------------------------------------
function x = bicgstab(theta,Q,PQ,QPQ,r,par)
% BiCGstab(ell) with preconditioning
% [x,rnrm] = bicgstab(theta,Q,PQ,QPQ,r,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=r 
% where Atilde=(I-Q*Q)*(A-theta*B)*(I-Q*Q').
% using (I-PQ*(QPQ\Q'))*inv(K) as preconditioner
%
% This function is specialized for use in JDQR.
% integer nmv: number of matrix multiplications
% rnrm: relative residual norm
%
%  par=[tol,mxmv,ell] where 
%    integer m: max number of iteration steps
%    real tol: residual reduction
%
% rnrm: obtained residual reduction
%
% -- References: ETNA

% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 1998, Gerard Sleijpen

% -- Initialization --
%

global Precond_Type

tol=par(1); max_it=par(2); l=par(3); n=size(r,1);
rnrm=1; nmv=0;
 
if max_it < 2 | tol>=1, x=r; return, end
%%% 0 step of bicgstab eq. 1 step of bicgstab
%%% Then x is a multiple of b

TP=Precond_Type;
if TP<1, 
if TP == 0, r=SkewProj(Q,PQ,QPQ,SolvePrecond(r)); end, tr=r;
else, tr=RepGS(Q,r); end
rnrm=norm(r); snrm=rnrm; tol=tol*snrm;

sigma=1; omega=1; 
x=zeros(n,1); u=zeros(n,1);
J1=2:l+1; 
   
%%% HIST=[0,1]; 

if TP <2 %% explicit preconditioning
% -- Iteration loop
while (nmv < max_it)

   sigma=-omega*sigma;
   for j = 1:l,
      rho=tr'*r(:,j);  bet=rho/sigma;
      u=r-bet*u;
      u(:,j+1)=mvp(theta,Q,PQ,QPQ,u(:,j));
      sigma=tr'*u(:,j+1);  alp=rho/sigma;
      r=r-alp*u(:,2:j+1);
      r(:,j+1)=mvp(theta,Q,PQ,QPQ,r(:,j)); nmv = nmv+2;
      x=x+alp*u(:,1);
      G(1,1)=r(:,1)'*r(:,1); rnrm=sqrt(G(1,1));
      if rnrm<tol, l=j; J1=2:l+1; r=r(:,1:l+1); break, end
   end
   
   for i=2:l+1 
     G(i,1:i)=r(:,i)'*r(:,1:i); G(1:i,i)=G(i,1:i)'; 
   end
   if TP, g=Q'*r; G=G-g'*g; end
   d=G(J1,1); gamma=G(J1,J1)\d;  
   rnrm=sqrt(real(G(1,1)-d'*gamma));   %%% compute norm in l-space
   %%%  HIST=[HIST;[nmv,rnrm/snrm]];

   x=x+r(:,1:l)*gamma;
   if rnrm < tol, break, end     %%% sufficient accuracy. No need to update r,u
   omega=gamma(l,1); gamma=[1;-gamma];
   u=u*gamma; r=r*gamma; 
   if TP, g=g*gamma; r=r-Q*g; end

   % rnrm = norm(r); 
end

else %% implicit preconditioning

I=eye(2*l); v0=I(:,1:l); s0=I(:,l+1:2*l);
y0=zeros(2*l,1); V=zeros(n,2*l);

while (nmv < max_it)

   sigma=-omega*sigma;
   y=y0; v=v0; s=s0;
   for j = 1:l,
      rho=tr'*r(:,j);  bet=rho/sigma;
      u=r-bet*u;
      if j>1,                 %%% collect the updates for x in l-space
         v(:,1:j-1)=s(:,1:j-1)-bet*v(:,1:j-1); 
      end
      [u(:,j+1),V(:,j)]=mvp(theta,Q,PQ,QPQ,u(:,j));
      sigma=tr'*u(:,j+1);  alp=rho/sigma;
      r=r-alp*u(:,2:j+1);
      if j>1, 
         s(:,1:j-1)=s(:,1:j-1)-alp*v(:,2:j); 
      end
      [r(:,j+1),V(:,l+j)]=mvp(theta,Q,PQ,QPQ,r(:,j));
      y=y+alp*v(:,1); 
      G(1,1)=r(:,1)'*r(:,1); rnrm=sqrt(G(1,1));
      if rnrm<tol, l=j; J1=2:l+1; s=s(:,1:l); break, end
   end
   nmv = nmv+2*l;

   for i=2:l+1 
     G(i,1:i)=r(:,i)'*r(:,1:i); G(1:i,i)=G(i,1:i)'; 
   end
   g=Q'*r; G=G-g'*g;                 %%% but, do the orth to Q implicitly
   d=G(J1,1); gamma=G(J1,J1)\d;   
   rnrm=sqrt(real(G(1,1)-d'*gamma)); %%% compute norm in l-space
   x=x+V*(y+s*gamma);

   %%% HIST=[HIST;[nmv,rnrm/snrm]];

   if rnrm < tol, break, end  %%% sufficient accuracy. No need to update r,u
   omega=gamma(l,1); gamma=[1;-gamma];
   u=u*gamma; r=r*gamma; 
   g=g*gamma; r=r-Q*g;        %%% Do the orth to Q explicitly
                              %%% In exact arithmetic not needed, but
                              %%% appears to be more stable.

end
end

if TP==1, x=SkewProj(Q,PQ,QPQ,SolvePrecond(x)); end
%%%  plot(HIST(:,1),log10(HIST(:,2)+eps),'*'), drawnow, pause
return

%----------------------------------------------------------------------
function [v,rnrm] = gmres0(theta,Q,PQ,QPQ,v,par)
% GMRES
% [x,rnrm] = gmres(theta,Q,Z,PQ,QPQ,v,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=b 
% where Atilde=(I-Q*Q)*(A-theta*B)*(I-Q*Q').
% using (I-PQ*(QPQ\Q'))*inv(K) as preconditioner
%
% If used as implicit preconditioner then FGMRES.
%
% par=[tol,m] where
%  integer m: degree of the minimal residual polynomial
%  real tol: residual reduction
%
% rnrm: obtained residual reduction
%
% -- References: Saad & Schultz SISC 1986

% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 1998, Gerard Sleijpen

% -- Initialization

global Precond_Type

tol=par(1); max_it=par(2); n = size(v,1);
rnrm = 1; j=0;

if max_it < 2 | tol>=1, return, end 
%%% 0 step of gmres eq. 1 step of gmres
%%% Then x is a multiple of b
 
H = zeros(max_it +1,max_it); Rot=[ones(1,max_it);zeros(1,max_it)];

TP=Precond_Type;

TP=Precond_Type; 
if TP<1
  if TP == 0,  v=SkewProj(Q,PQ,QPQ,SolvePrecond(v)); end
  rho0 = norm(v); v = v/rho0;
else
  v=RepGS(Q,v); 
end

V = [v];
tol = tol * rnrm; 
y = [ rnrm ; zeros(max_it,1) ];

while (j < max_it) & (rnrm > tol),
  j=j+1;
  [v,w]=mvp(theta,Q,PQ,QPQ,v); 
  if TP > 0
    if TP == 2, W=[W,w]; end 
    v=RepGS(Q,v,0); 
  end
  [v,h] = RepGS(V,v); H(1:size(h,1),j) = h;
  V = [V, v]; 
  for i = 1:j-1,
    a = Rot(:,i);
    H(i:i+1,j) = [a'; -a(2) a(1)]*H(i:i+1,j);
  end
  J=[j, j+1];
  a=H(J,j);
  if a(2) ~= 0
     cs = norm(a); 
     a = a/cs; Rot(:,j) = a;
     H(J,j) = [cs; 0];
     y(J) = [a'; -a(2) a(1)]*y(J);
  end 
  rnrm = abs(y(j+1));
end

J=[1:j];  
if TP == 2
  v = W(:,J)*(H(J,J)\y(J));
else
  v = V(:,J)*(H(J,J)\y(J));
end

if TP==1, v=SkewProj(Q,PQ,QPQ,SolvePrecond(v)); end

return
%%%======================================================================
function v = gmres(theta,Q,PQ,QPQ,v,par)
% GMRES 
% [x,nmv,rnrm] = gmres(theta,Q,PQ,QPQ,v,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=r 
% where Atilde=(I-Q*Q)*(A-theta*B)*(I-Q*Q').
% using (I-PQ*(QPQ\Q'))*inv(K) as preconditioner.
%
% If used as implicit preconditioner, then FGMRES.
%
% par=[tol,m] where
%  integer m: degree of the minimal residual polynomial
%  real tol: residual reduction
%
% nmv:  number of MV with Atilde
% rnrm: obtained residual reduction
%
% -- References: Saad
% Same as gmres0. However this variant uses MATLAB built-in functions
% slightly more efficient (see Sleijpen and van den Eshof).
%
%
% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 2002, Gerard Sleijpen

global Precond_Type

% -- Initialization
tol=par(1); max_it=par(2); n = size(v,1);
j=0;

if max_it < 2 | tol>=1, rnrm=1; return, end 
%%% 0 step of gmres eq. 1 step of gmres
%%% Then x is a multiple of b
 
H = zeros(max_it +1,max_it); Gamma=1; rho=1;

TP=Precond_Type; 
if TP<1
  if TP == 0,  v=SkewProj(Q,PQ,QPQ,SolvePrecond(v)); end
  rho0 = norm(v); v = v/rho0;
else
  v=RepGS(Q,v); rho0=1;
end

V = zeros(n,0); W=zeros(n,0);
tol0 = 1/(tol*tol); 
%% HIST=1; 

reduk=1; rho0=1; rho=1;

while (j < max_it) & (rho < tol0) 

  V=[V,v]; j=j+1;
  [v,w]=mvp(theta,Q,PQ,QPQ,v);
  if TP > 0
    if TP == 2, W=[W,w]; end 
    v=RepGS(Q,v,0); 
  end
  [v,h] = RepGS(V,v); 
  H(1:size(h,1),j)=h; gamma=H(j+1,j);
  
  if gamma==0, break %%% Lucky break-down
  else
    gamma= -Gamma*h(1:j)/gamma; 
    Gamma=[Gamma,gamma];
    rho=rho+gamma'*gamma;
  end     
          
red=rho0/rho; rho0=rho;
if (red>0.7+0.3*reduk),break,end
reduk=min(reduk,red);  

  %% HIST=[HIST;(gamma~=0)/sqrt(rho)]; 
    
end

if gamma==0; %%% Lucky break-down
   e1=zeros(j,1); e1(1)=rho0; rnrm=0; 
   if TP == 2
     v=W*(H(1:j,1:j)\e1); 
   else
     v=V*(H(1:j,1:j)\e1); 
   end 
else %%% solve in least square sense 
   e1=zeros(j+1,1); e1(1)=rho0; rnrm=1/sqrt(rho);
   if TP == 2
     v=W*(H(1:j+1,1:j)\e1); 
   else
     v=V*(H(1:j+1,1:j)\e1); 
   end 
end

if TP==1, v=SkewProj(Q,PQ,QPQ,SolvePrecond(v)); end
%% HIST=log10(HIST+eps); J=[0:size(HIST,1)-1]';
%% plot(J,HIST(:,1),'*'); drawnow,% pause
return

%======================================================================
%========== BASIC OPERATIONS ==========================================
%======================================================================
function v=MV(v)

  global Operator_A Operator_MVs

  if ischar(Operator_A)
    v = feval(Operator_A,v);
  else
    v = Operator_A*v;
  end

  Operator_MVs = Operator_MVs +size(v,2);
return
%----------------------------------------------------------------------
function [v,u]=mvp(theta,Q,Z,M,v)
% v=Atilde*v

global Precond_Type

  if Precond_Type > 0
    u=SkewProj(Q,Z,M,SolvePrecond(v)); 
    v=MV(u)-theta*u;
  else
    u=v; v=MV(u)-theta*u;
    v=SkewProj(Q,Z,M,SolvePrecond(v));
  end

return
%----------------------------------------------------------------------
function y=SolvePrecond(y)
% Action preconditioner

  global Precond_Form Precond_L Precond_U Precond_P Precond_Solves

    switch Precond_Form
      case 0, return
      case 1,     y=feval(Precond_L,y); 
      case 2,     y=feval(Precond_L,y,'preconditioner');
      case 3,     y=feval(Precond_L,y); 
                  y=feval(Precond_U,y);
      case 4,     y=feval(Precond_L,y,'L'); 
                  y=feval(Precond_L,y,'U');
      case 5,     y=Precond_L\y;
      case 6,     y=Precond_U\(Precond_L\y);
      case 7,     y=Precond_U\(Precond_L\(Precond_P*y));
    end
    Precond_Solves = Precond_Solves +size(y,2);
    
return
%----------------------------------------------------------------------
function  r=SkewProj(Q,Z,M,r);

   if ~isempty(Q), 
      r=r-Z*(M\(Q'*r));
   end 

return
%----------------------------------------------------------------------
function ppar=spar(par,nit)
% Changes par=[tol(:),max_it,ell]  to
% ppap=[TOL,max_it,ell] where 
% if lenght(tol)==1
%   TOL=tol
% else
%   red=tol(end)/told(end-1); tole=tol(end);
%   tol=[tol,red*tole,red^2*tole,red^3*tole,...]
%   TOL=tol(nit);
% end

k=size(par,2)-2;
ppar=par(1,k:k+2);

if k>1
   if nit>k
      ppar(1,1)=par(1,k)*((par(1,k)/par(1,k-1))^(nit-k));
   else
      ppar(1,1)=par(1,max(nit,1));
   end
end

ppar(1,1)=max(ppar(1,1),1.0e-8);

return
%
%======= Iterative methods for symmetric systems =======================
function x = cg(theta,Q,PQ,QPQ,r,par)
% CG
% [x,rnrm] = cg(theta,Q,PQ,QPQ,b,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=b 
% where Atilde=(I-PQ*QPQ^(-1)*Q')*(U\L\(A-theta)).
%
% par=[tol,m] where
%  integer m: degree of the minimal residual polynomial
%  real tol: residual reduction
%
% nmv:  number of MV with Atilde
% rnrm: obtained residual reduction
%
% -- References: Hestenes and Stiefel

% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 1998, Gerard Sleijpen

% -- Initialization

global Precond_Type

tol=par(1); max_it=par(2); n = size(r,1);
rnrm = 1; nmv=0; 

if max_it ==0 | tol>=1, x=r; return, end

x= zeros(n,1); u=zeros(n,1);

rho = norm(r); %% r=r/rho; rho=1; 
rho = rho*rho; snrm = rho; rho0=rho; reduk=1;
tol = tol*tol*rho; 
sigma=1;

%%%  HIST=rho;

while  ( rho > tol & nmv < max_it )

  if Precond_Type, 
    c=SkewProj(Q,PQ,QPQ,SolvePrecond(r)); rho=c'*r;
  else, c=r; end  

  beta=rho/sigma;   
  u=c-beta*u;
  c=MV(u) - theta*u; nmv=nmv+1;

  sigma=u'*c; alpha=rho/sigma;  
  x=x+alpha*u; 
  r=r-alpha*c; sigma=-rho;

  rho=r'*r;  

  %%%  HIST=[HIST;rho];

% if rho<snrm & rho<rho0, 
%   red=rho/rho0; reduk=min(reduk,red); 
%   if red>0.7+0.3*reduk, break, end
% end
% rho0=rho; 

end % while

%% HIST=log10(HIST+eps); J=[0:size(HIST,1)-1]';
%% plot(J,HIST(:,1),'*'); drawnow,% pause

return
%----------------------------------------------------------------------
function x = minres(theta,Q,Z,M,r,par)
% MINRES   with preconditioner K.
% [x,rnrm] = minres(theta,Q,Z,M,b,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=b 
% where Atilde=(I-Z*M^(-1)*Q')*(K\(A-theta)).
%
% A should be symmetric and K positive definit.
% Implicit preconditioning requires as many precondioner solves
% as explicit preconditioning.
%
% par=[tol,m] where
%  integer m: degree of the minimal residual polynomial
%  real tol: residual reduction
%
% nmv:  number of MV with Atilde
% rnrm: obtained residual reduction
%
% -- References: Paige and Saunders

% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 1998, Gerard Sleijpen

% -- Initialization

global Precond_Type
TP=Precond_Type;

tol=par(1); max_it=par(2); n = size(r,1);
rnrm = 1; nmv=0; 

if max_it ==0 | tol>=1, x=r; return, end

x=zeros(n,1); rho = norm(r); v = r/rho; snrm=rho;
beta = 0; v_old = zeros(n,1); 
beta_t = 0; c = -1; s = 0;
w = zeros(n,1);

tol=tol*rho;
%%% HIST = rho;

if TP < 0
  www = v; pv=v;
else
  pv=SkewProj(Q,Z,M,SolvePrecond(v)); 
  beta=sqrt(pv'*v); pv=pv/beta; v=v/beta; www=pv;
end

while  ( nmv < max_it  &  abs(rho) > tol )

   wv=RepGS(Q,MV(pv),0)-theta*pv; nmv=nmv+1;
   wv = wv - beta*v_old;
   alpha = pv'*wv; v_old = v; v = wv-alpha*v;
   if TP < 0
     beta = norm(v); v = v/beta; pv = v;
   else
     pv=SkewProj(Q,Z,M,SolvePrecond(v));
     beta=sqrt(pv'*v); pv=pv/beta; v=v/beta; 
   end

   l1 = s*alpha - c*beta_t; l2 = s*beta;

   alpha_t = -s*beta_t - c*alpha;  beta_t = c*beta;
   l0 = sqrt(alpha_t*alpha_t+beta*beta); 
   c = alpha_t/l0; s = beta/l0;

   ww = www - l1*w; www = pv - l2*w; w = ww/l0;

   x =  x + (rho*c)*w; rho =  s*rho; 

  %%%  HIST=[HIST;rho];
end % while

%% HIST=log10(HIST+eps); J=[0:size(HIST,1)-1]';
%% plot(J,HIST(:,1),'*'); drawnow, pause

return

%----------------------------------------------------------------------
function x = symmlq(theta,Q,Z,M,r,par)
% SYMMLQ
% [x,rnrm] = symmlq(theta,Q,Z,M,b,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=b 
% where Atilde=(I-Z*M^(-1)*Q')*(K\(A-theta)).
%
% A should be symmetric and K positive definit.
% Implicit preconditioning requires as many precondioner solves
% as explicit preconditioning..
%
% par=[tol,m] where
%  integer m: degree of the minimal residual polynomial
%  real tol: residual reduction
%
% nmv:  number of MV with Atilde
% rnrm: obtained residual reduction
%
% -- References: Paige and Saunders

% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 1998, Gerard Sleijpen

% -- Initialization

global Precond_Type
TP=Precond_Type;

tol=par(1); max_it=par(2); n = size(r,1);
rnrm = 1; nmv=0; 

if max_it ==0 | tol>=1, x=r; return, end

x=zeros(n,1); rho = norm(r); v = r/rho; snrm=rho;
beta = 0;  beta_t = 0; c = -1;  s = 0;
v_old = zeros(n,1);  gtt = rho; g = 0;
      
tol=tol*rho;
   
%%% HIST = rho;

if TP < 0
  pv=v; w=pv;
else
  pv=SkewProj(Q,Z,M,SolvePrecond(v)); 
  beta=sqrt(pv'*v); pv=pv/beta; v=v/beta; w=pv;
end

while  ( nmv < max_it  &  rho > tol )

   wv=RepGS(Q,MV(pv),0)-theta*pv; nmv=nmv+1;
   wv = wv - beta*v_old;
   alpha = pv'*wv; v_old = v; v = wv-alpha*v;
   if TP < 0
     beta = norm(v); v = v/beta; pv = v;
   else
     pv=SkewProj(Q,Z,M,SolvePrecond(v));
     beta=sqrt(pv'*v); pv=pv/beta; v=v/beta; 
   end

   l1 = s*alpha - c*beta_t; l2 = s*beta; 
      
   alpha_t = -s*beta_t - c*alpha; beta_t = c*beta;
   l0 = sqrt(alpha_t*alpha_t+beta*beta); 
   c = alpha_t/l0; s = beta/l0;

   gt = gtt - l1*g; gtt = -l2*g; g = gt/l0;

   rho = sqrt(gt*gt+gtt*gtt); 

   x = x + (g*c)*w + (g*s)*pv;
   w = s*w - c*pv; 
 
   %%%  HIST=[HIST;rho];
     
end % while

%% HIST=log10(HIST+eps); J=[0:size(HIST,1)-1]';
%% plot(J,HIST(:,1),'*'); drawnow, pause

return
%=======================================================================
%========== Orthogonalisation ==========================================
%=======================================================================
function [v,y]=RepGS(V,v,gamma)
% [v,y]=REP_GS(V,w)
% If V orthonormal then [V,v] orthonormal and w=[V,v]*y;
% If size(V,2)=size(V,1) then w=V*y;
%
% The orthonormalisation uses repeated Gram-Schmidt
% with the Daniel-Gragg-Kaufman-Stewart (DGKS) criterion.
%
% [v,y]=REP_GS(V,w,GAMMA)
% GAMMA=1 (default) same as [v,y]=REP_GS(V,w)
% GAMMA=0, V'*v=zeros(size(V,2)) and  w = V*y+v (v is not normalized).

 
% coded by Gerard Sleijpen, August 28, 1998

if nargin < 3, gamma=1; end

[n,d]=size(V);

if size(v,2)==0, y=zeros(d,0); return, end

nr_o=norm(v); nr=eps*nr_o; y=zeros(d,1);
if d==0
  if gamma, v=v/nr_o; y=nr_o; else, y=zeros(0,1); end, return
end

y=V'*v; v=v-V*y; nr_n=norm(v); ort=0;

while (nr_n<0.5*nr_o & nr_n > nr)
  s=V'*v; v=v-V*s; y=y+s; 
  nr_o=nr_n; nr_n=norm(v);     ort=ort+1; 
end

if nr_n <= nr, if ort>2, disp(' dependence! '), end
  if gamma  % and size allows, expand with a random vector
    if d<n, v=RepGS(V,rand(n,1)); y=[y;0]; else, v=zeros(n,0); end
  else, v=0*v; end
elseif gamma, v=v/nr_n; y=[y;nr_n]; end

return
%=======================================================================
%============== Sorts Schur form =======================================
%=======================================================================
function [Q,S]=SortSchur(A,sigma,gamma,kk)
%[Q,S]=SortSchur(A,sigma)
%  A*Q=Q*S with diag(S) in order prescribed by sigma.
%  If sigma is a scalar then with increasing distance from sigma.
%  If sigma is string then according to string
%  ('LM' with decreasing modulus, etc)
%
%[Q,S]=SortSchur(A,sigma,gamma,kk)
%  if gamma==0, sorts only for the leading element
%  else, sorts for the kk leading elements

  l=size(A,1);
  if l<2, Q=1;S=A; return, 
  elseif nargin==2, kk=l-1; 
  elseif gamma, kk=min(kk,l-1); 
  else, kk=1; sigma=sigma(1,:); end

%%%------ compute schur form -------------
  [Q,S]=schur(A); %% A*Q=Q*S, Q'*Q=eye(size(A));
%%% transform real schur form to complex schur form
  if norm(tril(S,-1),1)>0, [Q,S]=rsf2csf(Q,S); end

%%%------ find order eigenvalues ---------------
  I = SortEig(diag(S),sigma); 

%%%------ reorder schur form ----------------
  [Q,S] = SwapSchur(Q,S,I(1:kk)); 

return
%----------------------------------------------------------------------
function I=SortEig(t,sigma);
%I=SortEig(T,SIGMA) sorts the indices of T.
%
% T is a vector of scalars, 
% SIGMA is a string or a vector of scalars.
% I is a permutation of (1:LENGTH(T))' such that:
%   if SIGMA is a vector of scalars then
%   for K=1,2,...,LENGTH(T) with KK = MIN(K,SIZE(SIGMA,1))
%      ABS( T(I(K))-SIGMA(KK) ) <= ABS( T(I(J))-SIGMA(KK) ) 
%      SIGMA(kk)=INF: ABS( T(I(K)) ) >= ABS( T(I(J)) ) 
%         for all J >= K

if ischar(sigma)
  switch sigma
    case 'LM'
      [s,I]=sort(-abs(t));
    case 'SM'
      [s,I]=sort(abs(t));
    case 'LR';
      [s,I]=sort(-real(t));
    case 'SR';
      [s,I]=sort(real(t));
    case 'BE';
      [s,I]=sort(real(t)); I=twistdim(I,1);
  end
else

  [s,I]=sort(abs(t-sigma(1,1))); 
  ll=min(size(sigma,1),size(t,1)-1);
  for j=2:ll
    if sigma(j,1)==inf
      [s,J]=sort(abs(t(I(j:end)))); J=flipdim(J,1);
    else
      [s,J]=sort(abs(t(I(j:end))-sigma(j,1)));
    end
    I=[I(1:j-1);I(J+j-1)];
  end 

end

return
%----------------------------------------------------------------------
function t=twistdim(t,k)

  d=size(t,k); J=1:d; J0=zeros(1,2*d);
  J0(1,2*J)=J; J0(1,2*J-1)=flipdim(J,2); I=J0(1,J);
  if k==1, t=t(I,:); else, t=t(:,I); end

return
%----------------------------------------------------------------------
function [Q,S]=SwapSchur(Q,S,I)
% [Q,S]=SwapSchur(QQ,SS,P)
%    QQ and SS are square matrices of size K by K
%    P is the first part of a permutation of (1:K)'.
%
%    If    M = QQ*SS*QQ'  and  QQ'*QQ = EYE(K), SS upper triangular
%    then  M*Q = Q*S      with   Q'*Q = EYE(K),  S upper triangular
%    and   D(1:LENGTH(P))=DD(P) where D=diag(S), DD=diag(SS)
%
%    Computations uses Givens rotations.

  kk=min(length(I),size(S,1)-1);
  j=1; while (j<=kk & j==I(j)), j=j+1; end; 
  while j<=kk
    i=I(j);
    for k=i-1:-1:j
      q = [S(k,k)-S(k+1,k+1),S(k,k+1)];
      if q(1) ~= 0
        q = q/norm(q);
        G = [[q(2);-q(1)],q'];
        J = [k,k+1];
        Q(:,J) = Q(:,J)*G;
        S(:,J) = S(:,J)*G;
        S(J,:) = G'*S(J,:);
      end
      S(k+1,k) = 0;
    end
    I=I+(I<i);   
    j=j+1; while (j<=kk & j==I(j)), j=j+1; end
  end

return
%----------------------------------------------------------------------
function [Q,Z,S,T]=SortQZ(A,B,sigma,gamma,kk)
%
% [Q,Z,S,T]=SORTQZ(A,B,SIGMA)
%   A and B are K by K matrices, SIGMA is a complex scalar or string.
%   SORTQZ computes the qz-decomposition of (A,B) with prescribed
%   ordering: A*Q=Z*S, B*Q=Z*T; 
%             Q and Z are K by K unitary,
%             S and T are K by K upper triangular.
%   The ordering is as follows:
%   (DAIG(S),DAIG(T)) are the eigenpairs of (A,B) ordered
%   as prescribed by SIGMA.
%

% coded by Gerard Sleijpen, version Januari 12, 1998


  l=size(A,1); 
  if l<2; Q=1; Z=1; S=A; T=B; return
  elseif nargin==3, kk=l-1; 
  elseif gamma, kk=min(kk,l-1); 
  else, kk=1; sigma=sigma(1,:); end

%%%------ compute qz form ----------------
  [S,T,Z,Q]=qz(A,B); Z=Z'; S=triu(S);
  
%%%------ sort eigenvalues ---------------

  I=SortEigPair(diag(S),diag(T),sigma); 
 
%%%------ sort qz form -------------------
  [Q,Z,S,T]=SwapQZ(Q,Z,S,T,I(1:kk)); 

return
%----------------------------------------------------------------------
function I=SortEigPair(s,t,sigma)
% I=SortEigPair(S,T,SIGMA)
%   S is a complex K-vectors, T a positive real K-vector
%   SIGMA is a string or a vector of pairs of complex scalars.
%   SortEigPair gives the index set I that sorts the pairs (S,T).
%
%   If SIGMA is a pair of scalars then the sorting is 
%   with increasing "chordal distance" w.r.t. SIGMA.
%
%  The chordal distance D between a pair A and a pair B is defined as follows.
%  Scale A by a scalar F such that NORM(F*A)=1 and F*A(2)>=0,
%  scale B by a scalar G such that NORM(G*B)=1 and G*B(2)>=0,
%  then D(A,B)=ABS((F*A)*RROT(G*B)) where RROT(alpha,beta)=(beta,-alpha)


% coded by Gerard Sleijpen, version Januari 14, 1998


n=sign(t); n=n+(n==0); t=abs(t./n); s=s./n; 

if ischar(sigma)
  switch sigma
    case {'LM','SM'}
    case {'LR','SR','BE'}
      s=real(s);
  end
  [s,I]=sort((-t./sqrt(s.*conj(s)+t.*t)));
  switch sigma
    case {'LM','LR'}
      I=flipdim(I,1);
    case {'SM','SR'}
    case 'BE'
      I=twistdim(I,1);  
  end
else

  n=sqrt(sigma.*conj(sigma)+1); ll=size(sigma,1); 
  tau=[ones(ll,1)./n,-sigma./n]; tau=tau.';

  n=sqrt(s.*conj(s)+t.*t); s=[s./n,t./n];

  [t,I]=sort(abs(s*tau(:,1))); 
  ll = min(ll,size(I,1)-1); 
  for j=2:ll
    [t,J]=sort(abs(s(I(j:end),:)*tau(:,j))); 
    I=[I(1:j-1);I(J+j-1)];
  end 

end

return

%----------------------------------------------------------------------
function I=SortEigPairVar(S,T,gamma)
% Compute harmonic Ritz vectors y; and Rayleigh quotients theta
% Depending on gamma, sort on harmonic Ritz value, Rayleigh quotient,
% residual norm or ...

  [X,L]=eig(S,T); L=S*X; X=T*X;
  M=diag(X'*X);
  R=diag(L'*L); % |product theta and htheta|
  theta=M.*R;   % square of theta 
  Res=R-theta;  % square of residual norm of (y,theta)
  htheta=M.\R;  % square of harm Ritz value

  if gamma==0, [d,I]=sort(htheta); return, end
  if gamma==1, [d,I]=sort(theta); return, end
  if gamma==2, [d,I]=sort(Res); return, end
  if gamma==3, [d,I]=sort(R); return, end

return


%----------------------------------------------------------------------
function [Q,Z,S,T]=SwapQZ(Q,Z,S,T,I)
% [Q,Z,S,T]=SwapQZ(QQ,ZZ,SS,TT,P)
%    QQ and ZZ are K by K unitary,  SS and TT are K by K uper triangular.
%    P is the first part of a permutation of (1:K)'.
%
%    Then Q and Z are K by K unitary, S and T are K by K upper triangular,
%    such that, for A = ZZ*SS*QQ' and B = ZZ*T*QQ', we have 
%    A*Q = Z*S, B*Q = Z*T  and LAMBDA(1:LENGTH(P))=LLAMBDA(P) where 
%    LAMBDA=DIAG(S)./DIAGg(T) and LLAMBDA=DIAG(SS)./DIAG(TT).
%
%    Computation uses Givens rotations. 
%

% coded by Gerard Sleijpen, version October 12, 1998
  

  kk=min(length(I),size(S,1)-1);
  j=1; while (j<=kk & j==I(j)), j=j+1; end
  while j<=kk
    i=I(j);
    for k = i-1:-1:j, 
      %%% i>j, move ith eigenvalue to position j 
      J = [k,k+1]; 
      q = T(k+1,k+1)*S(k,J) - S(k+1,k+1)*T(k,J);
      if q(1) ~= 0 
        q = q/norm(q);
        G = [[q(2);-q(1)],q'];
        Q(:,J) = Q(:,J)*G; 
        S(:,J) = S(:,J)*G; T(:,J) = T(:,J)*G;
      end 
      if abs(S(k+1,k))<abs(T(k+1,k)), q=T(J,k); else q=S(J,k); end
      if q(2) ~= 0
        q=q/norm(q);
        G = [q';q(2),-q(1)];
        Z(:,J) = Z(:,J)*G'; 
        S(J,:) = G*S(J,:); T(J,:) = G*T(J,:);
      end 
      T(k+1,k) = 0;
      S(k+1,k) = 0; 
    end
    I=I+(I<i); 
    j=j+1; while (j<=kk & j==I(j)), j=j+1; end
  end

return

%=======================================================================
%======= SET PARAMETERS ================================================
%=======================================================================
function [n,nselect,sigma,SCHUR,...
         jmin,jmax,tol,maxit,V,INTERIOR,SHOW,PAIRS,JDV,FIX,OLD,...
         lsolver,par] = ReadOptions(varargin)
% Read options and set defaults

global Operator_A Precond_Form Precond_L Precond_U Precond_P Precond_Type

Operator_A = varargin{1};  

%%% determine dimension
if ischar(Operator_A)
  n=-1;
  if exist(Operator_A) ~=2
    msg=sprintf('  Can not find the M-file ''%s.m''  ',Operator_A);
    errordlg(msg,'MATRIX'),n=-2;
  end
  if n==-1, eval('n=feval(Operator_A,[],''dimension'');','n=-1;'), end
else
  [n,n] = size(Operator_A);
  if any(size(Operator_A) ~= n)
    msg=sprintf('  The operator must be a square matrix or a string.  ');
    errordlg(msg,'MATRIX'),n=-3;
  end
end

%%% defaults
nselect0= 5; 
SCHUR   = 0;
jmin    = -1;
jmax    = -1;
p0      = 5; % jmin=nselect+p0
p1      = 5; % jmax=jmin+p1
tol     = 1e-8; 
maxit   = 200;
V       = zeros(0,0);
INTERIOR= 0;
SHOW    = 0;
PAIRS   = 0;
JDV     = 0;
OLD     = 0.0001;
FIX     = 1000;
lsolver = 'gmres';
ls_maxit= 200; 
ls_tol  = [1,0.7];  
ell     = 4;
par     = [ls_tol,ls_maxit,ell];
TP      = 0;
V0      = 'ones(n,1)+0.1*rand(n,1)';  %%%% 'v'

options=[]; sigma=[]; varg=[]; Precond_L = []; Precond_U = []; Precond_P = [];
for j = 2:nargin
  if isstruct(varargin{j})
    options = varargin{j};
  elseif ischar(varargin{j}) 
    if length(varargin{j}) == 2 & isempty(sigma)
      sigma = varargin{j};
    elseif isempty(Precond_L)
      Precond_L=varargin{j};
    elseif isempty(Precond_U)
      Precond_U=varargin{j};
    end
  elseif length(varargin{j}) == 1
    varg = [varg,varargin{j}];
  elseif min(size(varargin{j}))==1 
    sigma = varargin{j}; if size(sigma,1)==1, sigma=conj(sigma'); end 
  elseif isempty(Precond_L)
    Precond_L=varargin{j};
  elseif isempty(Precond_U)
    Precond_U=varargin{j};
  elseif isempty(Precond_P)
    Precond_P=varargin{j};
  end
end

%------- set targets ----------------
if ischar(sigma)
  sigma0=sigma; sigma=upper(sigma);
  switch sigma
    case {'LM','LR','SR','BE','SM'}
    otherwise, sigma=ShiftOrPrecond(sigma0);
  end
end

[s,I]=sort(varg); I=flipdim(I,2);  
J=[]; j=0; 
while j<length(varg)
  j=j+1; jj=I(j); s=varg(jj);
  if isreal(s) & (s == fix(s)) & (s > 0)
    if n==-1
      n=s; eval('v=feval(Operator_A,zeros(n,0));','n=-1;')
      if n>-1, J=[J,jj]; end
    end
  else
    if ischar(sigma),  sigma=ShiftOrPrecond(sigma); end
    if isempty(sigma), sigma=s; end 
    J=[J,jj];  
  end 
end
varg(J)=[];

if n==-1,
  msg1=sprintf('  Cannot find the dimension of ''%s''.  \n',Operator_A);
  msg2=sprintf('  Put the dimension n in the parameter list:  \n  like');
  msg3=sprintf('\t\n\n\t jdqr(''%s'',n,..),  \n\n',Operator_A);
  msg4=sprintf('  or let\n\n\t n = %s(',Operator_A);
  msg5=sprintf('[],''dimension'')\n\n  give n.');
  msg=[msg1,msg2,msg3,msg4,msg5];
  errordlg(msg,'MATRIX')
end

nselect=[]; 
if n<2, return, end

if length(varg) == 1
   nselect=min(n,varg);
elseif length(varg)>1
   if isempty(sigma), sigma=varg(end); varg(end)=[]; end
   nselect=min(n,min(varg));
end


%------- number of eigs to be computed ----------
if isempty(nselect), nselect=min(n,nselect0); end

fopts = []; if ~isempty(options), fopts=fieldnames(options); end

%------- preconditioner -------------------------
Precond_L=findfield(options,fopts,'pr',[]);
[L,ok]=findfield(options,fopts,'l_',Precond_L);
if ok & ~isempty(Precond_L),
   msg =sprintf('A preconditioner is defined in');
   msg =[msg,sprintf('\n''Precond'', but also in ''L_precond''.')];
   msg=[msg,sprintf('\nWhat is the correct one?')];
   button=questdlg(msg,'Preconditioner','L_Precond','Precond','L_Precond');
   if strcmp(button,'L_Precond'), 
     Precond_L = L;
   end
else
   Precond_L = L;
end

if ~isempty(Precond_L)
  Precond_U=findfield(options,fopts,'u_',[]);
  Precond_P=findfield(options,fopts,'p_',[]);
end

if isempty(Precond_L), ls_tol  = [0.7,0.49]; end

TP=findfield(options,fopts,'ty',TP);
n=SetPrecond(n,TP); if n<0, return, end

%------- max, min dimension search subspace ------
jmin=min(n,findfield(options,fopts,'jmi',jmin));
jmax=min(n,findfield(options,fopts,'jma',jmax));
if jmax < 0
   if jmin<0, jmin=min(n,nselect+p0); end
   jmax=min(n,jmin+p1); 
else
   if jmin<0, jmin=max(1,jmax-p1); end
end 
maxit=abs(findfield(options,fopts,'ma',maxit));

%------- initial search subspace ----------------
V=findfield(options,fopts,'v',[]);
[m,d]=size(V); 
if m~=n 
  if m>n, V = V(1:n,:); end
  if m<n, V = [V;0.001*rand(n-m-1,d)]; end
end
V=orth(V); [m,d]=size(V);
if d==0, nr=0; while nr==0, V=eval(V0); nr=norm(V); V=V/nr; end, end

INTERIOR=findfield(options,fopts,'te',INTERIOR,[0,1,1.1,1.2],{'st','ha'});

[sigma,n]=ShiftIfInterior(sigma,INTERIOR,n); if n<0, return, end

%------- Other options --------------------------
tol   = findfield(options,fopts,'to',tol);
SCHUR = findfield(options,fopts,'sch',SCHUR,[0,1,0.5]);
PAIRS = findfield(options,fopts,'pai',PAIRS,[0,1]);
JDV   = findfield(options,fopts,'av',JDV,[0,1]);
SHOW  = findfield(options,fopts,'di',SHOW,[0,1,2]);
OLD   = max(abs(findfield(options,fopts,'tr',OLD,[0,OLD,inf])),10*tol);
FIX   = max(abs(findfield(options,fopts,'fix',0,[0,FIX,inf])),0);
lsolver= lower(findfield(options,fopts,'lso',lsolver));

%------- linear solver --------------------------
switch lsolver
  case {'exact'}
    Precond_L=[]; Precond_Form=0;
    if ischar(Operator_A)
  msg=sprintf('The operator must be a matrix for ''exact''.');
  msg=[msg,sprintf('\nDo you want to solve the correction equation')];
  msg=[msg,sprintf('\naccurately with an iterative solver (BiCGstab)?')];
      button=questdlg(msg,'Solving exactly','Yes','No','Yes');
      if strcmp(button,'No'), n=-1; return, end
    end
  case {'iluexact'}
    if ischar(Precond_L)
  msg=sprintf('The preconditioner must be matrices for ''iluexact''.');
      errordlg(msg,'Solving with ''iluexact''')
      n=-1; return
    end
  case {'olsen'}
  case {'cg'}
    ls_tol=1.0e-12;
    if Precond_Form>0, Precond_Type=2; end
  case {'minres','symmlq'}
    ls_tol=1.0e-12;
    if Precond_Form>0, Precond_Type=2; end
  case 'gmres'
    ls_maxit=5;
  case {'bicgstab','cgstab'}
    lsolver='bicgstab'; ls_tol=1.0e-10;
  otherwise
    error(['Unknown method ''' lsolver '''.']);
end

ls_maxit= abs(findfield(options,fopts,'ls_m',ls_maxit));
ls_tol  = abs(findfield(options,fopts,'ls_t',ls_tol));
ell     = round(findfield(options,fopts,'ls_e',ell));
par=[ls_tol,ls_maxit,ell];

if SHOW
   ShowChoices(n,nselect,sigma,SCHUR,jmin,jmax,tol,maxit,V,...
         INTERIOR,SHOW,PAIRS,JDV,FIX,OLD,lsolver,ls_tol,ls_maxit,ell)
end

return
%-------------------------------------------------------------------
function sigma=ShiftOrPrecond(sigma);

global Precond_L Precond_U

 if isempty(sigma) | ~ischar(sigma), return, end
 if ~isempty(Precond_L) & ~ischar(Precond_L), return, end
 if exist(sigma)==2 & (isempty(Precond_L) | isempty(Precond_U))
   ok=1; eval('v=feval(sigma,zeros(n,1));','ok=0')
   if ok
     Precond_U=Precond_L; Precond_L=sigma; sigma=[]; 
   end
 end

return
%-------------------------------------------------------------------
function n=SetPrecond(n,TP)
% finds out how the preconditioners are defined (Precond_Form)
% and checks consistency of the definitions.
%
% If M is the preconditioner then P*M=L*U. Defaults: L=U=P=I.
%
% Precond_Form
%       0:   no L
%       1:   L M-file, no U,     L ~= A
%       2:   L M-file, no U,     L == A
%       3:   L M-file, U M-file, L ~= A, U ~= A, L ~=U
%       4:   L M-file, U M-file, L == U
%       5:   L matrix, no U
%       6:   L matrix, U matrix  no P
%       7:   L matrix, U matrix, P matrix
%
% Precond_Type
%       0:   Explicit left
%       1:   Explicit right
%       2:   Implicit
%

  global Operator_A ...
         Precond_Form Precond_Solves ...
         Precond_Type ...
         Precond_L Precond_U Precond_P

  Precond_Type = 0;
  if ischar(TP)
    TP=lower(TP(1,1));
    switch TP
      case 'l'
        Precond_Type = 0;
      case 'r'
        Precond_Type = 1;
      case 'i'
        Precond_Type = 2;   
    end
  else
    Precond_Type=max(0,min(fix(TP),2));
  end

  Precond_Solves = 0;

  % Set type preconditioner
  Precond_Form=0;
  if isempty(Precond_L), Precond_Type=-1; return, end

  if ~isempty(Precond_U) & ischar(Precond_L)~=ischar(Precond_U)
    msg=sprintf('  L and U should both be strings or matrices');
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end
  if ~isempty(Precond_P) & (ischar(Precond_P) | ischar(Precond_L))
    msg=sprintf('  P can be specified only if P, L and U are matrices'); 
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end  
  tp=1+4*~ischar(Precond_L)+2*~isempty(Precond_U)+~isempty(Precond_P);
  if tp==1, tp = tp + strcmp(Precond_L,Operator_A); end
  if tp==3, tp = tp + strcmp(Precond_L,Precond_U); end
  if tp==3 & strcmp(Precond_U,Operator_A)
    msg1=sprintf('  If L and A use the same M-file,')
    msg2=sprintf('\n  then so should U.'); 
    errordlg([msg1,msg2],'PRECONDITIONER'), n=-1; return
  end
  if tp>5, tp=tp-1; end, Precond_Form=tp;

  % Check consistency definitions
  if tp<5 & exist(Precond_L) ~=2
    msg=sprintf('  Can not find the M-file ''%s.m''  ',Precond_L); 
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end

  ok=1; 
  if tp == 2
    eval('v=feval(Operator_A,zeros(n,1),''preconditioner'');','ok=0;')
    if ~ok
       msg='Preconditioner and matrix use the same M-file';
       msg1=sprintf(' %s.   \n',Precond_L);
       msg2='Therefore the preconditioner is called';
       msg3=sprintf(' as\n\n\tw=%s(v,''preconditioner'')\n\n',Precond_L); 
       msg4='Put this "switch" in the M-file.';
       msg=[msg,msg1,msg2,msg3,msg4];
    end
  end

  if tp == 4 | ~ok
    ok1=1;
    eval('v=feval(Precond_L,zeros(n,1),''L'');','ok1=0;')
    eval('v=feval(Precond_L,zeros(n,1),''U'');','ok1=0;')
    if ok1 
      Precond_Form = 4; Precond_U = Precond_L; ok=1;
    else
      if tp == 4
        msg='L and U use the same M-file';
        msg1=sprintf(' %s.m   \n',Precond_L);
        msg2='Therefore L and U are called';
        msg3=sprintf(' as\n\n\tw=%s(v,''L'')',Precond_L); 
        msg4=sprintf(' \n\tw=%s(v,''U'')\n\n',Precond_L); 
        msg5=sprintf('Check the dimensions and/or\n');
        msg6=sprintf('put this "switch" in %s.m.',Precond_L);
        msg=[msg,msg1,msg2,msg3,msg4,msg5,msg6]; 
      end
      errordlg(msg,'PRECONDITIONER'), n=-1; return
    end
  end

  if tp==1 | tp==3
    eval('v=feval(Precond_L,zeros(n,1));','ok=0')
    if ~ok
       msg=sprintf('''%s'' should produce %i-vectors',Precond_L,n); 
       errordlg(msg,'PRECONDITIONER'), n=-1; return 
    end
  end

  if tp==3
    if exist(Precond_U) ~=2
      msg=sprintf('  Can not find the M-file ''%s.m''  ',Precond_U);
      errordlg(msg,'PRECONDITIONER'), n=-1; return
    else
      eval('v=feval(Precond_U,zeros(n,1));','ok=0')
      if ~ok
        msg=sprintf('''%s'' should produce %i-vectors',Precond_U,n);
        errordlg(msg,'PRECONDITIONER'), n=-1; return 
      end
    end
  end
  
  if tp==5
    if min([n,2*n]==size(Precond_L)) 
      Precond_U=Precond_L(:,n+1:2*n); Precond_L=Precond_L(:,1:n); 
      Precond_Form=6;
    elseif min([n,3*n]==size(Precond_L)) 
      Precond_U=Precond_L(:,n+1:2*n); Precond_P=Precond_L(:,2*n+1:3*n);
      Precond_L=Precond_L(:,1:n); Precond_Form=7;
    elseif ~min([n,n]==size(Precond_L)) 
      msg=sprintf('The preconditioning matrix\n');
      msg2=sprintf('should be %iX%i or %ix%i ([L,U])\n',n,n,n,2*n); 
      msg3=sprintf('or %ix%i ([L,U,P])\n',n,3*n); 
      errordlg([msg,msg2,msg3],'PRECONDITIONER'), n=-1; return
    end
  end

  if tp==6 & ~min([n,n]==size(Precond_L) & [n,n]==size(Precond_U))
    msg=sprintf('Both L and U should be %iX%i.',n,n); n=-1;
    errordlg(msg,'PRECONDITIONER'), n=-1; return 
  end

  if tp==7 & ~min([n,n]==size(Precond_L) & ...
       [n,n]==size(Precond_U) & [n,n]==size(Precond_P))
    msg=sprintf('L, U, and P should all be %iX%i.',n,n); n=-1;
    errordlg(msg,'PRECONDITIONER'), n=-1; return 
  end

return
%-------------------------------------------------------------------
function [sigma,n]=ShiftIfInterior(sigma,INTERIOR,n)

  if isempty(sigma)
    if INTERIOR, sigma=0; else, sigma = 'LM'; end, return
  elseif ischar(sigma)
    switch sigma
      case {'LM','LR','SR','BE'}
        if INTERIOR == 0, return, end
      case {'SM'}
        sigma=0; return
      otherwise
        if INTERIOR, sigma=0; else, sigma='LM'; end, return
     end
  else, return
  end

   msg1=sprintf('\n   The choice sigma = ''%s'' does not match',sigma);
   msg2=sprintf('\n   the search for INTERIOR eigenvalues.');
   msg3=sprintf('\n   Specify a numerical value for sigma');
   msg4=sprintf(',\n   for instance, a value that is ');
   switch sigma
      case {'LM'}
         % sigma='SM';
         msg5=sprintf('absolute large.\n');
         msg4=[msg4,msg5];
      case {'LR'}
         % sigma='SP'; % smallest positive real
         msg5=sprintf('positive large.\n');
         msg4=[msg4,msg5];
      case {'SR'}
         % sigma='LN'; % largest negative real
         msg5=sprintf('negative and absolute large.\n');
         msg4=[msg4,msg5];
      case {'BE'}
         % sigma='AS'; % alternating smallest pos., largest neg
         msg4=sprintf('???.\n');
   end
   msg=[msg1,msg2,msg3,msg4];
   msg5=sprintf('   Do you want to continue with sigma=0?');
   msg=[msg,msg5];
   button=questdlg(msg,'Finding Interior Eigenvalues','Yes','No','Yes');
   if strcmp(button,'Yes'), sigma=0, else, n=-1; end

return
%------------------------------------------------------------------------
function x = boolean(x,gamma,string)
%Y = BOOLEAN(X,GAMMA,STRING)
%  GAMMA(1) is the default. 
%  If GAMMA is not specified, GAMMA = 0.
%  STRING is a cell of accepted strings. 
%  If STRING is not specified STRING = {'n' 'y'}
%  STRING{I} and GAMMA(I) are accepted expressions for X 
%  If X=GAMMA(I) then Y=X. If the first L characters
%  of X matches those of STRING{I}, then Y=GAMMA(I+1).
%  Here L=SIZE({STRING{1},2).
%  For other values of X, Y=GAMMA(1);

if nargin < 2, gamma=[0,0,1]; end
if nargin < 3, string={'n' 'y'}; end

if ischar(x)
  l=size(string{1},2);
  i=min(find(strncmpi(x,string,l)));
  if isempty(i), i=0; end, x=gamma(i+1);
elseif max((gamma-x)==0)
elseif gamma(end) == inf
else, x=gamma(1);
end
  
return
%------------------------------------------------------------------------
function [a,ok]=findfield(options,fopts,str,default,gamma,stri)
% Searches the fieldnames in FOPTS for the string STR.
% The field is detected if only the first part of the fieldname
% matches the string STR. The search is case insensitive.
% If the field is detected, then OK=1 and A is the fieldvalue.
% Otherwise OK=0 and A=DEFAULT

   l=size(str,2); j=min(find(strncmpi(str,fopts,l)));
   if ~isempty(j)
      a=getfield(options,char(fopts(j,:))); ok=1;
      if nargin == 5, a = boolean(a,[default,gamma]); 
      elseif nargin == 6, a = boolean(a,[default,gamma],stri); end
   elseif nargin>3
      a=default; ok=0;
   else
      a=[]; ok=0;
   end

return
%===========================================================================
%============= OUTPUT FUNCTIONS ============================================
%===========================================================================
function ShowOptions
 fprintf('\n')
 fprintf('PROBLEM\n')
 fprintf('            A: [ square matrix | string ]\n');
 fprintf('      nselect: [ positive integer {5} ]\n\n');
 fprintf('TARGET\n')
 fprintf('        sigma: [ scalar | row or vector of scalars |\n');
 fprintf('                 ''LM'' | {''SM''} | ''LR'' | ''SR'' | ''BE'' ]\n\n');

 fprintf('OPTIONS\n');
 fprintf('        Schur: [ yes | {no} ]\n');
 fprintf('          Tol: [ positive scalar {1e-8} ]\n');
 fprintf('         Disp: [ yes | {no} | 2 ]\n');
 fprintf('         jmin: [ positive integer {nselect+5} ]\n');
 fprintf('         jmax: [ positive integer {jmin+5} ]\n');
 fprintf('        MaxIt: [ positive integer {200} ]\n');
 fprintf('           v0: [ size(A,1) by p vector of scalars {rand(size(A,1),1)} ]\n');
 fprintf('    TestSpace: [ Standard | {Harmonic} ]\n');
 fprintf('        Pairs: [ yes | {no} ]\n');
 fprintf('    AvoidStag: [ yes | {no} ]\n');
 fprintf('        Track: [ {yes} | no | non-negative scalar {1e-4} ]\n');
 fprintf('     FixShift: [ yes | {no} | non-negative scalar {1e+3} ]\n');
 fprintf('      LSolver: [ {gmres} | bicgstab ]\n');
 fprintf('       LS_Tol: [ row of positive scalars {[1,0.7]} ]\n');
 fprintf('     LS_MaxIt: [ positive integer {5} ]\n');
 fprintf('       LS_ell: [ positive integer {4} ]\n');
 fprintf('      Precond: ');
 fprintf('[ n by n matrix {identity} | n by 2n matrix | string ]\n');
 fprintf('    L_Precond: same as ''Precond''\n');
 fprintf('    U_Precond: [ n by n matrix {identity} | string ]\n');
 fprintf('    P_Precond: [ n by n matrix {identity} ]\n');
 fprintf(' Type_Precond: [ {left} | right | implicit ]\n');    
if 1
 fprintf('\nPRECONDITIONER as input argument\n');
 fprintf('            M: [ n by n or n by 2n matrix {identity} | string ]\n');
 fprintf('or M=L*U (strings or matrices), or M=P\\(L*U) (matrices) where\n');
 fprintf('            L: [ n by n matrix {identity} | string ]\n');
 fprintf('            U: [ n by n matrix {identity} | string ]\n');
 fprintf('            P: [ n by n {identity} ]\n');
end
fprintf('\n')
return
%------------------------------------------------------------------------
function ShowChoices(n,nselect,sigma,SCHUR,...
         jmin,jmax,tol,maxit,V,INTERIOR,SHOW,PAIRS,JDV,FIX,OLD,...
         lsolver,ls_tol,ls_maxit,ell)

global Operator_A

  fprintf('\n'),fprintf('PROBLEM\n')
  fprintf('%13s: %s\n','A',StrOp(Operator_A));
  fprintf('%13s: %i\n','dimension',n);
  fprintf('%13s: %i\n','nselect',nselect);

  fprintf('\mTARGET\n')
  if ischar(sigma)
    fprintf('%13s: ''%s''\n','sigma',sigma)
  else 
    Str=ShowLambda(sigma);
    fprintf('%13s: %s\n','sigma',Str)
  end

  fprintf('\nOPTIONS\n');
  fprintf('%13s: %g\n','Schur',SCHUR)
  fprintf('%13s: %g\n','Tol',tol)
  fprintf('%13s: %i\n','Disp',SHOW)
  fprintf('%13s: %i\n','jmin',jmin)
  fprintf('%13s: %i\n','jmax',jmax)
  fprintf('%13s: %i\n','MaxIt',maxit)
  fprintf('%13s: %s\n','v0',StrOp(V))
 %fprintf('%13s: %g\n','TestSpace',INTERIOR);
  fprintf('%13s: %i\n','Pairs',PAIRS)
  fprintf('%13s: %i\n','AvoidStag',JDV)
  fprintf('%13s: %g\n','Track',OLD)
  fprintf('%13s: %g\n','FixShift',FIX)
  fprintf('%13s: ''%s''\n','LSolver',lsolver)


  switch lsolver
    case {'exact','iluexact','olsen'}
    case {'cg','minres','symmlq','gmres','bicgstab'}
    if length(ls_tol)>1
      str=sprintf('%g ',ls_tol);
      fprintf('%13s: [ %s]\n','LS_Tol',str)
    else
      fprintf('%13s: %g\n','LS_Tol',ls_tol)
    end
    fprintf('%13s: %i\n','LS_MaxIt',ls_maxit)
    if strcmp(lsolver,'bicgstab')
      fprintf('%13s: %i\n','LS_ell',ell)
    end
  end

  DisplayPreconditioner(n);

  string1='\n%13s: ''%s'''; string2='\n\t       %s';
  switch INTERIOR
  case 0
      fprintf(string1,'TestSpace','Standard, W = V ')
     % fprintf(string2,'V W: V orthogonal')
     % fprintf(string2,'W=A*V')
   case 1
      fprintf(string1,'TestSpace','Harmonic, W = A*V - sigma*V')
     % fprintf(string2,'V W: V and W orthogonal')
     % fprintf(string2,'AV-Q*E=W*R where AV=A*V-sigma*V and E=Q''*AV')
   case 1.1
      fprintf(string1,'TestSpace','Harmonic, W = A*V - sigma*V')
      fprintf(string2,'V W AV: V and W orthogonal')
      fprintf(string2,'AV=A*V-sigma*V, AV-Q*Q''*AV=W*R')
   case 1.2
      fprintf(string1,'TestSpace','Harmonic, W = A*V - sigma*V')
      fprintf(string2,'V W: W orthogonal')
      fprintf(string2,'W=AV-Q*E where AV=A*V-sigma*V and E=Q''*AV')
  otherwise
      fprintf(string1,'TestSpace','Experimental')
  end % switch INTERIOR
  fprintf('\n\n')

return
%------------------------------------------------------------------------
function msg=StrOp(Op)

  if ischar(Op)
    msg=sprintf('''%s''',Op);
  elseif issparse(Op), [n,k]=size(Op);
    msg=sprintf('[%ix%i sparse]',n,k);
  else, [n,k]=size(Op);
    msg=sprintf('[%ix%i double]',n,k);
  end

return
%------------------------------------------------------------------------
function DisplayPreconditioner(n)
  global Precond_Form Precond_Type ...
         Precond_L Precond_U Precond_P

  FP=Precond_Form;
  switch Precond_Form
    case 0,    
      fprintf('%13s: %s\n','Precond','No preconditioner'); 
    case {1,2,5}
      fprintf('%13s: %s','Precond',StrOp(Precond_L));
      if FP==2, fprintf('  (M\\v = %s(v,''preconditioner''))',Precond_L); end 
      fprintf('\n')
    case {3,4,6,7} 
      fprintf('%13s: %s\n','L precond',StrOp(Precond_L));
      fprintf('%13s: %s\n','U precond',StrOp(Precond_U));
      if FP==7, fprintf('%13s: %s\n','P precond',StrOp(Precond_P)); end
      if FP==4,fprintf('%15s(M\\v = %s(%s(v,''L''),''U''))\n',...
          '',Precond_L,Precond_L); end 
  end

  if FP
  switch Precond_Type
    case 0, str='explicit left';
    case 1, str='explicit right';
    case 2, str='implicit';
  end
  fprintf('%15sTo be used as %s preconditioner.\n','',str)
  end

return
%-------------------------------------------------------------------
function varargout=ShowLambda(lambda,kk)

for k=1:size(lambda,1);
  if k>1, Str=[Str,sprintf('\n%15s','')]; else, Str=[]; end
  rlambda=real(lambda(k,1)); ilambda=imag(lambda(k,1));
  Str=[Str,sprintf(' %+11.4e',rlambda)]; 
  if abs(ilambda)>100*eps*abs(rlambda) 
    if ilambda>0                    
      Str=[Str,sprintf(' + %10.4ei',ilambda)];
    else
      Str=[Str,sprintf(' - %10.4ei',-ilambda)];
    end
  end
end

if nargout == 0
  if nargin == 2
    Str=[sprintf('\nlambda(%i) =',kk),Str];
  else
    Str=[sprintf('\nDetected eigenvalues:\n\n%15s',''),Str];
  end
  fprintf('%s\n',Str)
else
  varargout{1}=Str;
end

return
%===========================================================================
function DispResult(s,nr,gamma)

  if nargin<3, gamma=0; end

  extra='';
  if nr > 100*eps & gamma
     extra='    norm > 100*eps !!! ';
  end

  if gamma<2 | nr>100*eps
     fprintf('\n %35s: %0.5g\t%s',s,nr,extra)
  end

return
%===========================================================================
function  STATUS0=MovieTheta(n,nit,Lambda,jmin,tau,LOCKED,SHRINK)
% MovieTheta(n,nit,Lambda,jmin,tau,nr<t_tol,j==jmax);

global Operator_A Rschur EigMATLAB CIRCLE MovieAxis M_STATUS

f=256;

if nargin==0, 
  if ~isempty(CIRCLE)
    figure(f), buttons(f,-2); hold off, refresh
  end
  return 
end

if nit==0
  EigMATLAB=[];
  if ~ischar(Operator_A) & n<201
    EigMATLAB=eig(full(Operator_A));
  end
  CIRCLE=0:0.005:1; CIRCLE=exp(CIRCLE*2*pi*sqrt(-1)); 
  if ischar(tau)
    switch tau
    case {'LR','SR'}
      if ~ischar(Operator_A)
        CIRCLE=norm(Operator_A,'inf')*[sqrt(-1),-sqrt(-1)]+1;
      end
    end
  end

  Explanation(~isempty(EigMATLAB),f-1)
end

if gcf~=f, figure(f), end
jmin=min(jmin,length(Lambda));
plot(real(Lambda(1:jmin)),imag(Lambda(1:jmin)),'bo'); 
if nit>0, axis(MovieAxis), end, hold on
pls='bo'; if SHRINK, pls='mo'; end
plot(real(Lambda(jmin+1:end)),imag(Lambda(jmin+1:end)),pls)
plot(real(EigMATLAB),imag(EigMATLAB),'cp'); 
THETA=diag(Rschur); plot(real(THETA),imag(THETA),'k*');

x=real(Lambda(1)); y=imag(Lambda(1));
plot(x,y,'kd'), pls='ks';
if LOCKED, plot(x,y,'ks'), pls='ms'; end

if ischar(tau), tau=0; end
delta=Lambda([1,jmin])-tau; 
if length(CIRCLE)~=2, delta=abs(delta); 
  plot(real(tau),imag(tau),pls)
else,  delta=real(delta); end
for i=1:1+SHRINK,  
  zeta=delta(i)*CIRCLE+tau;
  plot(real(zeta),imag(zeta),'r:'), 
end

if nit==0
  buttons(f,-1); hold off, zoom on
end
title('legend see figure(255)')
if SHRINK, STATUS0=buttons(f,2); if STATUS0, return, end, end
STATUS0=buttons(f,1); drawnow, hold off

return

%=============================================================
function Explanation(E,f)

  if gcf~=f, figure(f), end

  HL=plot(0,0,'kh',0,0,'bo'); hold on
  StrL=str2mat('Detected eigenvalues','Approximate eigenvalues');
  if E
    HL=[HL;plot(0,0,'cp')];
    StrL=str2mat(StrL,'Exact eigenvenvalues');
  end
  HL=[HL;plot(0,0,'kd',0,0,'ks',0,0,'r:')];
 % StrL=str2mat(StrL,'Tracked app. eig.','target',...
 %      'inner/outer bounds for restart');
 StrL=str2mat(StrL,'Selected approximate eigenvalue','target',...
       'inner/outer bounds for restart');

  legend(HL,StrL), hold on, drawnow, hold off
  title('legend for figure(256)')

return

%=============================================================
function STATUS0=buttons(f,push)
% push=0: do nothing
% push>0: check status buttons,
%                 STATUS0=1 if break else STATUS0=0.
% push=-1: make buttons for pause and break
% push=-2: remove buttons

global M_STATUS MovieAxis

  if push>0 % check status buttons

    ud = get(f,'UserData'); 
    if ud.pause ==1, M_STATUS=1; end 
    if ud.break ==1, STATUS0=1;
      ud.pause=0; ud.break=0; set(f,'UserData',ud); return, 
    else, STATUS0=0; end   

   if push>1
     ud.pause=0; set(f,'UserData',ud); 
     while M_STATUS
        ud = get(f,'UserData'); pause(0.1) 
        if ud.pause, M_STATUS=0; end
        if ud.break, M_STATUS=0; STATUS0=1; end
        MovieAxis=axis; 
     end
     ud.pause=0; ud.break=0; set(f,'UserData',ud);
   end


  elseif push==0, STATUS0=0; return
  elseif push==-1 % make buttons

    ud = [];
    h = findobj(f,'Tag','pause');
    if isempty(h)
      ud.pause = 0;
      pos = get(0,'DefaultUicontrolPosition');
      pos(1) = pos(1) - 15;
      pos(2) = pos(2) - 15;
      str = 'ud=get(gcf,''UserData''); ud.pause=1; set(gcf,''UserData'',ud);';
      uicontrol( ...
          'Style','push', ...
          'String','Pause', ...
          'Position',pos, ...
          'Callback',str, ...
          'Tag','pause');
    else
      set(h,'Visible','on');              % make sure it's visible
      if ishold
        oud = get(f,'UserData');
        ud.pause = oud.pause;             % don't change old ud.pause status
      else
        ud.pause = 0;
      end
    end
    h = findobj(f,'Tag','break');
    if isempty(h)
      ud.break = 0;
      pos = get(0,'DefaultUicontrolPosition');
      pos(1) = pos(1) + 50;
      pos(2) = pos(2) - 15;
      str = 'ud=get(gcf,''UserData''); ud.break=1; set(gcf,''UserData'',ud);';
      uicontrol( ...
          'Style','push', ...
          'String','Break', ...
          'Position',pos, ...
          'Callback',str, ...
          'Tag','break');
    else
      set(h,'Visible','on');            % make sure it's visible
      if ishold
        oud = get(f,'UserData');
        ud.break = oud.break;           % don't change old ud.break status
      else
        ud.break = 0;
      end
    end
    set(f,'UserData',ud); M_STATUS=0;

    STATUS0=0; hold off, zoom on

    MA=axis; nl=MA/10;
    MovieAxis = MA-max(nl([2,4])-nl([1,3]))*[1,-1,1,-1];

  else % remove buttons

    set(findobj(f,'Tag','pause'),'Visible','off');
    set(findobj(f,'Tag','break'),'Visible','off');
    STATUS0=0; return, refresh

  end
return
