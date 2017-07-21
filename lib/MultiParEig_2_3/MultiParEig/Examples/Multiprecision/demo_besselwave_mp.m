%DEMO_BESSELWAVE_MP  Bessel wave system is solved as a two-parameter eigenvalue problem 
% with Chebyshev collocation and multiprecision refinement
% 
% This example reconstructs Table 1  
% Lew Yan Voon & Willatzen, Helmholtz equation in parabolic rotational 
% coordinates: application to wave problems in quantum mechanics and
% acoustics, Math. Comp. Sim. 65 (2004) 337--349.
%
% See also: DEMO_BESSELWAVE1, DEMO_BESSELWAVE2, DEMO_BESSELWAVE3, DEMO_BESSELWAVE_FIGS, BESSELWAVE_MEP, TWOPAREIGS

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 2.1.2017

% Example output:
% 
% First 20 modes for Dirichlet & Dirichlet b.c. in multiprecision
% ===============================================================
% p  |    lambda                |    omega  
% ----------------------------------------------------------
% 0  | -0.0000000000000000000000 |  4.8096511153915455372433 
% 1  | -0.0000000000000000000000 |  6.2831853071795864769253 
% 2  | -0.0000000000000000000000 |  7.6634119404150246312289 
% 0  | 13.4667958177820646479784 |  7.8727664021316195496864 
% 3  | -0.0000000000000000000000 |  8.9868189158181283506158 
% 1  | 21.7319156456917796315528 |  9.3564714089201262312410 
% 4  | -0.0000000000000000000000 | 10.2712446036813651126028 
% 2  | 29.6901295549229639287235 | 10.7728606336823131748283 
% 0  | 39.9742137058085642349441 | 10.8920989630208841825531 
% 0  | -0.0000000000000000000000 | 11.0401562205726212991932 
% 5  | -0.0000000000000000000000 | 11.5269183937890995828129 
% 3  | 37.6483315799820638489124 | 12.1419230182830022314445 
% 1  | 56.2696972787430322763073 | 12.3619967607514745452369 
% 1  | -0.0000000000000000000000 | 12.5663706143591729538506 
% 6  | -0.0000000000000000000000 | 12.7603237918479670124732 
% 4  | 45.6939627970971804428420 | 13.4755195579225769057820 
% 2  | 72.1395035833624179050866 | 13.7835409995494130797138 
% 0  | 81.0097881199184116837157 | 13.8703996099951480164153 
% 7  | -0.0000000000000000000000 | 13.9758640010010399180307 
% 2  |  0.0000000000000000000000 | 14.0311733396312375070741 
% ----------------------------------------------------------

is_numeric_type_supported('mp'); % check if MCT is installed

N = 60;
XiSq = 1;
EtaSq = 1;

% We compute enough modes (8 for p=0,..,8) and then select the smallest sqrt(mu)
% We select only modes with lambda>0 (since some are close to zero and
% negative, we take lambda>-1)

% -------------------------------------------------------------------
% Table 1 - Dirichlet & Dirichlet
tic
modes = [];
for p=0:8
    fprintf('Computing modes for p=%d in double precision\n',p);
    [A1,B1,C1,A2,B2,C2,z1,z2,G1,G2] = besselwave_mep(N,N,p,sqrt(XiSq),sqrt(EtaSq),[1 1]);
    [lambda,mu,X1,X2] = twopareigs(A1,B1,C1,A2,B2,C2,8);
    modes = [modes; p*ones(length(lambda),1) lambda sqrt(mu)];
end

posmodes = modes(find(real(modes(:,2))>-1),:);
[tilda,ord] = sort(real(posmodes(:,3)));
t1 = toc

Table1 = posmodes(ord(1:20),:);  % first 20 modes
fprintf('\nFirst 20 modes for Dirichlet & Dirichlet b.c. in double precision\n');
disp('===================================================================');
disp('p  |    lambda    |    omega  ');
disp('------------------------------');
for k=1:20
    fprintf('%1i  | %11.8f | %11.8f \n',posmodes(ord(k),1),posmodes(ord(k),2),posmodes(ord(k),3))
end
fprintf('------------------------------\n\n');

% -------------------------------------------------------------------
% Table 1 - Dirichlet & Dirichlet in quadruple precision
tic 
modes2 = mp([]);
opts = [];
opts.fp_type = 'mp';
opts.double_init = 1;
opts.refine_mp = 10;
for p=0:8
    fprintf('Computing modes for p=%d with TRQI refinement in multiprecision\n',p);
    [A1,B1,C1,A2,B2,C2,z1,z2,G1,G2] = besselwave_mep(N,N,p,sqrt(XiSq),sqrt(EtaSq),[1 1],opts);
    [lambda2,mu2] = twopareigs(A1,B1,C1,A2,B2,C2,8,opts);
    modes2 = [modes2; p*ones(length(lambda2),1) lambda2 sqrt(mu2)];
end
t2 = toc

posmodes2 = modes2(find(real(modes2(:,2))>-1),:);
[tilda,ord2] = sort(real(posmodes2(:,3)));

Table1 = posmodes2(ord2(1:20),:);  % first 20 modes
fprintf('\nFirst 20 modes for Dirichlet & Dirichlet b.c. in multiprecision\n');
disp('===============================================================');
disp('p  |    lambda                |    omega  ');
disp('----------------------------------------------------------');
for k=1:20
    fprintf('%1i  | %25.22f | %25.22f \n',posmodes2(ord2(k),1),posmodes2(ord2(k),2),posmodes2(ord2(k),3))
end
fprintf('----------------------------------------------------------\n\n');

