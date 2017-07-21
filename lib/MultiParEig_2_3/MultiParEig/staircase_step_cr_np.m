function [Delta, full, rep] = staircase_step_cr_np(Delta, opts)

%STAIRCASE_STEP_CR_NP  column-row reduction step of staircase algorithm for Delta matrices
% 
% [Delta, full, rep] = STAIRCASE_STEP_CR_NP(Delta, opts)
% reduces Delta{1} to full column rank using column and row compression on 
% (Delta{1}, Delta{2}), ... ,(Delta{1}, Delta{k+1}), where k = length(Delta)-1
%
% Options in opts:
%   - rank1 : -1 (if different from -1 this is fixed rank for 1st reduction)
%   - rank2 : -1 (if different from -1 this is fixed rank for 2nd reduction)
%   - rrqr (0) : set to 1 to use rank-revealing qr with column pivoting instead of SVD
%
% Output:
%  - Delta : reduceed Delta matrices
%  - full : final matrix is full rank (1) or no (0)
%  - rep : data with details on the rank selection
%
% See also: STAIRCASE_STEP_RC_NP, EXTRACT_REGULAR_PART_NP.

% Reference: A. Muhic, B. Plestenjak: On the quadratic two-parameter 
% eigenvalue problem and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 16.09.2016 : rank1 and rank2 options for fixed rank choices, report
% BP 23.09.2016 : option of using qr with column pivoting instead of SVD
% PH 26.11.2016 : support for arbitrary numeric types

% Last revision: 26.11.2016

if nargin<2, opts=[]; end
if isfield(opts,'rank1'), rank1 = opts.rank1; else, rank1 = -1; end
if isfield(opts,'rank2'), rank2 = opts.rank2; else, rank2 = -1; end
if isfield(opts,'rrqr'),  rrqr = opts.rrqr;   else, rrqr = 0;   end
 
k = length(Delta)-1; % number of parameters
class_t = superiorfloat(Delta{:});

opts.call = 'RedColRow 1';
opts.firstlastgap = 2;
opts.fixedrank = rank1;
opts.rrqr = rrqr*2;
[ran1, tilda1, tilda2, V1, extra1] = extended_svd(Delta{1}, opts); %#ok<*ASGLU>
rep = [1 size(Delta{1}) extra1];
% if the selection of the fixed rank failed, we do not fix the second rank
if rep(end)>=0
    rank2 = -1;
end

% column nullity
if ran1 == size(Delta{1},2)
    full = 1;
    return;
else
    full = 0;
end

POM = numeric_t([],class_t);
for i = 2:k+1
    Delta{i} = Delta{i}*V1;
    POM = [POM Delta{i}(:,ran1+1:end)]; %#ok<*AGROW>
end

opts.call = 'RedColRow 2';
opts.firstlastgap = 1;
opts.fixedrank = rank2;
opts.rrqr = rrqr*1;
[ran2, U2, tilda1, tilda2, extra2] = extended_svd(POM, opts);
rep = [rep; 2 size(POM) extra2];

if ran2 == size(POM,1)
    % reduction leads to empty matrices
    for i = 1:k+1
        Delta{i} = numeric_t([],class_t); 
    end
    return
end

Delta{1} = U2(:,ran2+1:end)'*Delta{1}*V1(:,1:ran1);
for i = 2:k+1
    Delta{i} = U2(:,ran2+1:end)'*Delta{i}(:,1:ran1);
end

