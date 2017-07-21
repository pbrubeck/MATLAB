% DEMO_LAME_EXTENDED  Lame system is solved as a two-parameter eigenvalue problem 
% with Chebyshev collocation
%
% This example reconstructs Table 1 and Table 2 (page 246) from 
% Morrison & Lewis, Charge singularity at the corner of a flat plate,
% SIAM J. Appl. Math. 31 (1976) 233--250
% and demonstrates convergence with respect to number of collocation points
%
% See also: DEMO_LAME, LAME_MEP, TWOPAREIGS, DEMO_LAME_CONV

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

format long

t = [1:7 9:15]/8; % values [0.125 ... 1.875] in table 1
Table1 = t(:);

for N=40:20:100
    value = [];
    for k=1:length(t)
        [A1,B1,C1,A2,B2,C2] = lame_mep(t(k)*pi,N,N);
        A1=A1+5*B1; A2=A2+5*B2; % we shift because A2 is singular
        [lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,1);
        value(k) = max(roots([1 1 -mu])); % positive root of z^2+z-mu=0
    end
    Table1 = [Table1 value(:)];
end

disp('====================================================================================================');
disp('Table 1 in Morrison and Lewis, SIAM J. Appl. Math. 31 (1976) 233--250');
disp('----------------------------------------------------------------------------------------------------');
disp('   t                   N=40                N=60                N=80                N=100');
disp('----------------------------------------------------------------------------------------------------');
disp(Table1)

format short e
disp('----------------------------------------------------------------------------------------------------');
disp('Differences to the N=100');
disp('----------------------------------------------------');
disp('   t            N=40         N=60         N=80 ');
disp('----------------------------------------------------');
disp([Table1(:,1) Table1(:,2)-Table1(:,5) Table1(:,3)-Table1(:,5) Table1(:,4)-Table1(:,5)])

t = [0.04021 0.11610 0.16670 0.22432 0.28858 0.875 0.900 0.925 0.950 1.050 1.075 1.100 1.125 1.875 1.900 1.925 1.950];
Table2 = t(:);

for N=40:20:100
    value = [];
    for k=1:length(t)
        [A1,B1,C1,A2,B2,C2] = lame_mep(t(k)*pi,N,N);
        A1=A1+5*B1; A2=A2+5*B2; % we shift because A2 is singular
        [lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,1);
        value(k) = max(roots([1 1 -mu])); % positive root of z^2+z-mu=0
    end
    Table2 = [Table2 value(:)];
end

format long
disp('====================================================================================================');
disp('Table 2 in Morrison and Lewis, SIAM J. Appl. Math. 31 (1976) 233--250');
disp('----------------------------------------------------------------------------------------------------');
disp('   t                   N=40                N=60                N=80                N=100');
disp('----------------------------------------------------------------------------------------------------');
disp(Table2)

format short e
disp('----------------------------------------------------------------------------------------------------');
disp('Differences to the N=100');
disp('----------------------------------------------------');
disp('   t            N=40         N=60         N=80 ');
disp('----------------------------------------------------');
disp([Table2(:,1) Table2(:,2)-Table2(:,5) Table2(:,3)-Table2(:,5) Table2(:,4)-Table2(:,5)])
