function [Delta, full] = staircase_step_cr_np(Delta, opts)

%STAIRCASE_STEP_CR_NP  column-row reduction step of staircase algorithm for Delta matrices
% 
% [Delta, full] = STAIRCASE_STEP_CR_NP(Delta, opts)
% reduces Delta{1} to full column rank using column and row compression on 
% (Delta{1}, Delta{2}), ... ,(Delta{1}, Delta{k+1}), where 
% k = length(Delta)-1
%
% See also: STAIRCASE_STEP_RC_NP, EXTRACT_REGULAR_PART_NP.

% Reference: A. Muhic, B. Plestenjak: On the quadratic two-parameter 
% eigenvalue problem and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<2, opts=[]; end
k = length(Delta)-1; % number of parameters

opts.call = 'RedColRow 1';
opts.firstlastgap = 2;
[ran1, tilda1, tilda2, V1] = extended_svd(Delta{1}, opts);

% column nullity
if ran1 == size(Delta{1},2)
    full = 1;
    return;
else
    full = 0;
end;

POM = [];
for i = 2:k+1
    Delta{i} = Delta{i}*V1;
    POM = [POM Delta{i}(:,ran1+1:end)];
end

opts.call = 'RedColRow 2';
opts.firstlastgap = 1;
[ran2, U2, tilda1, tilda2] = extended_svd(POM, opts);
if ran2 == size(POM,1)
    % reduction leads to empty matrices
    for i = 1:k+1
        Delta{i} = []; 
    end
    return
end

Delta{1} = U2(:,ran2+1:end)'*Delta{1}*V1(:,1:ran1);
for i = 2:k+1
    Delta{i} = U2(:,ran2+1:end)'*Delta{i}(:,1:ran1);
end

