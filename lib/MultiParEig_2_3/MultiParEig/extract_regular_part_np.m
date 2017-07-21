function [Delta,report] = extract_regular_part_np(Delta,opts)

%EXTRACT_REGULAR_PART_NP  Regular projections of Delta matrices from singular MEP 
%
% [Delta,report] = EXTRACT_REGULAR_PART_NP(Delta,opts)
% returns "Delta" matrices of the projected regular problem
% with D{i} = P*Delta{i}*Q, i = 1, ..., k+1 such that D{1} is of full rank, 
% matrices (D{1})^(-1)*D{2}, ..., (D{1})^(-1)*D{k+2} commute (hopefully)
% and P and Q are maximal possible matrices with orthogonal columns. It
% uses staircase algorithm with column-row and row-column compressions.
%
% opts:
%  - ranksequence : [] (m x 3 matrix with rows [i rank1 rank2], where 
%    i = 1 for cr and 2 for rc reduction, rank1 and rank2 are fixed ranks)
%  - mingap : (1) minimum gap (log10) between singular values at rank
%
% Output:
%  - Delta : reduced "Delta" matrices
%  - report with details on the rank selection and the staircase algorithm

% Reference: A. Muhic, B. Plestenjak: On the quadratic two-parameter 
% eigenvalue problem and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 16.09.2016 : ranksequence option added, oneway option removed, report

% Last revision: 16.9.2016

if nargin<2, opts=[]; end
if isfield(opts,'ranksequence'), RS = opts.ranksequence; else RS = [];      end
if isfield(opts,'mingap'),       mingap = opts.mingap;   else mingap = 1;   end

report = [];

for k = 1:size(RS,1)
    opts.rank1 = RS(k,2);
    opts.rank2 = RS(k,3);
    if RS(k,1) == 1
        [Delta, full, rep] = staircase_step_cr_np(Delta, opts); %#ok<*ASGLU>
    else
        [Delta, full, rep] = staircase_step_rc_np(Delta, opts);
    end
    report = [report; rep];
    if log10(min(rep(:,6)./rep(:,7)))<mingap
        return
    end
    opts = rmfield(opts,'rank1');
    opts = rmfield(opts,'rank2');
    if report(end,end) >=0 
        break
    end
end

[n1,n2] = size(Delta{1});
if ~isempty(report) 
    lastrank = report(end,4);
    if report(end,end)>=0 % continue with the staircase, not rank sequence
        lastrank = 0;
    end
else
    lastrank = 0;
end

% if we end with the square pencil in the first phase, this is the end
if (n1 ~= n2) || (lastrank==0)
    
    full = 0;
    while ~full
        [Delta, full, rep] = staircase_step_cr_np(Delta, opts);
        report = [report; rep]; %#ok<*AGROW>
        if log10(min(rep(:,6)./rep(:,7)))<mingap
            return
        end
    end
    
    % if we have a square pencil with nonsingular Delta{1} after the column-row compressions, we are done
    [n1,n2] = size(Delta{1});
    if ~isempty(report)
        lastrank = report(end,4);
    else
        lastrank = 0;
    end
    full = (n1 == n2) && (n1 == lastrank);
    % if not we continue with the row-column compression phase
    while ~full
        [Delta, full, rep] = staircase_step_rc_np(Delta, opts);
        report = [report; rep];
        if log10(min(rep(:,6)./rep(:,7)))<mingap
            return
        end
    end
    
end