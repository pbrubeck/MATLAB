%DEMO_MATHIEU_JD_MP   demo for method twopareigs_jd in multiprecision
%
% This example computes first eigenmodes for the Mathieu two-parameter eigenvalue problem
% using Jacobi-Davidson method for two-parameter eigenvalue problems in
% multiprecision
%
% See also: DEMO_MATHIEU_JD

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 4.12.2016

n = 100;   % size of the matrices
neig = 3; % number of wanted eigenvalues

%% First version is Jacobi-Davidson (double precision)
[A1,B1,C1,A2,B2,C2] = mathieu_mep(n,n,2,2,1); 
opts = [];
opts.target = [0 0];
opts.M1 = inv(A1-opts.target(1)*B1-opts.target(2)*C1);
opts.M2 = inv(A2-opts.target(1)*B2-opts.target(2)*C2);
opts.harmonic = 0;
opts.extraction = 'minmu';
opts.innersteps = 0;
opts.minsize = 6;
opts.maxsize = 12;
opts.maxmaxsize = 20;
opts.maxsteps = 500;
opts.delta = 10*eps*max(norm(A1),norm(A2));
opts.window = 0;
opts.showinfo = 2;
opts.forcereal = 1;
opts.refine = 4;
opts
tic; [lambda1,mu1,X1,Y1,X1l,Y1l,conv,hist] = twopareigs_jd(A1,B1,C1,A2,B2,C2,neig,opts); t1 = toc
mu1

%% Second version is Jacobi-Davidson (multi precision)

opts = [];
opts.fp_type = 'mp';
[A1,B1,C1,A2,B2,C2] = mathieu_mep(n,n,2,2,1,opts); 
opts = [];
opts.target = [0 0];
opts.M1 = inv(A1-opts.target(1)*B1-opts.target(2)*C1);
opts.M2 = inv(A2-opts.target(1)*B2-opts.target(2)*C2);
opts.harmonic = 0;
opts.extraction = 'minmu';
opts.innersteps = 0;
opts.minsize = 6;
opts.maxsize = 12;
opts.maxmaxsize = 20;
opts.maxsteps = 500;
opts.delta = 1e3*mp('eps')*max(norm(A1),norm(A2));
opts.switcheps = 1e15*opts.delta; % note that switcheps/delta is much larger as in double precision, here we switch to TRQI earlier
opts.window = 0;
opts.showinfo = 1;
opts.forcereal = 1;
opts.refine = 5;
opts
tic; [lambda2,mu2,X2,Y2,X2l,Y2l,conv,hist] = twopareigs_jd(A1,B1,C1,A2,B2,C2,neig,opts); t2 = toc
mu2
