%DEMO_TWOPAREIGS_JD   demo for method twopareigs_jd

% In this numerical example we generate two-parameter eigenvalue problems 
% with sparse 400 x 400 matrices and we are looking for the eigenvalues closest 
% to the target [0 0].  
%
% For numerical solutions we use Jacobi-Davidson method in twopareigs_jd
%
% See also: TWOPAREIGS_JD, RANDOM_SPARSE_2EP.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 21.11.2016

% Preparation of input data 
 
% We generate a MEP problem is such way that the solutions are known and we can
% compare the obtained results with the true results.
n = 400;
[A1,B1,C1,A2,B2,C2,lambda,mu] = random_sparse_2ep(n,0.03); 
target = [0 0];  		 % we are looking for the closest eigenvalue to the target

M1 = inv(full(A1-target(1)*B1-target(2)*C1)); % preconditioners
M2 = inv(full(A2-target(1)*B2-target(2)*C2)); 


% Setting up parameters

solveM1 = @(x) M1*x;  % we can pass preconditioners as function handles
solveM2 = @(x) M2*x;

tic
OPTS=[];							
OPTS.M1 = solveM1;    % preconditioners 
OPTS.M2 = solveM2;    
OPTS.minsize = 5;     
OPTS.maxsize = 10;
OPTS.maxsteps = 1000;   % maximum number of outer iterations
OPTS.extraction = 'mindist'; 
OPTS.reschange = 0; % 10^(-6); % we change to minimal residual when residual is less then epschange
OPTS.innersteps = 2;   % number of GMRES steps
OPTS.innertol = 1e-15;  % tolerance for the GMRES method
OPTS.target = target;
OPTS.delta = 5e-6;
OPTS.showinfo = 2; % set to 1 to see more information
OPTS.harmonic = 1;
OPTS.window = 0;

neig = 20;

% -----------------------------------------------------------------------
% First computation returns 10 eigenvalues
[ll,ee,Xr,Yr,Xl,Yl,res] = twopareigs_jd(A1,B1,C1,A2,B2,C2,10,OPTS);
toc

% -----------------------------------------------------------------------
% We use the computed eigenvectors as input for the second computation
% and compute 25 new eigenvalues (selection test prevents the method to
% choose already computed eigenvector again) 
% Although this works, we do not get so good eigenvalues in the second
% phase. It is better to compute all eigenvalues in one run. The reason for
% this is that the subspaces in the second run do not contain vectors close
% to the eigenvectors we are looking for.
OPTS.XPr = Xr;
OPTS.YPr = Yr;
OPTS.XPl = Xl;
OPTS.YPl = Yl;

OPTS.window = 10;
OPTS.maxsize = 20;
OPTS.maxsteps = 500;

tic
[ll2,ee2,Xr2,Yr2,Xl2,Yl2,res2] = twopareigs_jd(A1,B1,C1,A2,B2,C2,neig-5,OPTS);
toc

% we check the indices of the computed eigenvalues to see 
% if the computed eigenvalues are indeed closest to the target 

% first we sort exact eigenvalues on their distance from the target
dist = abs(lambda-target(1)).^2+abs(mu-target(2)).^2;
[tmp,ord] = sort(dist);
lambda = lambda(ord);
mu = mu(ord);

n1 = length(ll);
ind1 = [];
for k=1:n1
    dist = abs(lambda-ll(k)).^2+abs(mu-ee(k)).^2;
    [tmp,pos] = min(dist);
    ind1(k) = pos;
end
if n1==0
    ind1 = 0;
end

n2 = length(ll2);
ind2 = [];
for k=1:n2
    dist = abs(lambda-ll2(k)).^2+abs(mu-ee2(k)).^2;
    [tmp,pos] = min(dist);
    ind2(k) = pos;
end
if n2==0
    ind2 = 0;
end

ind1 % indices of the first group
ind2 % indices of the second group



