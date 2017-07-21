function [r,choice]= numrank(d, opts)

%NUMRANK   Find a gap in singular values
%
% [r,choice] = NUMRANK(d, opts) returns rank from a set of singular values d
%
% Arguments: 
%   d : a set of singular values (e.g., result of the function svd)
%   opts: options
% 
% Options in opts:
%   - rankeps : 1e-12 (discard all singular values smaller from rankeps)
%   - heuristic : 1 (take the largest gap if it is clear, set to 0 to take the number singular values above d(1)*rankeps
%   - mingapdif : 0.5 (minimal difference between the largest and the second largest gap, where gap(i) = log10(d(i)/d(i+1))
%   - firstlastgap : 0 (special strategy for selecting rank - default is 0)
%       1 : if gap is not clear, select first gap - use when it is safer to detect smaller rank and you know that the rank can not be full
%       2 : if gap is not clear, select last gap - use when it is safer to detect larger rank
%   - fixedrank : -1 (if different from -1 then select fixedrank if gap(fixedrank) is larger than fixedrankgap) 
%   - fixedrankgap : 1 (required minimal gap at fixedrank)
%   - showrank : 0 (display diagnostics, 1: partial, 2: detail)
%
% Output:
%   r : numerical rank
%   choice : how was rank selected (0: larger than rankeps*d(1), -1: fixedrank, 1: fixedrank overruled, 2,3,4,5,6: heuristic version)

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 16.09.2015: new options fixedrank and fixedrankgap, report choice, mingap -> mingapdif
% PH 22.11.2016: precision-independent version
% BP 26.11.2016: minor change in the printout
% PH 26.11.2016: code simplifications and clean-ups.

% Last revision: 26.11.2016
class_t = class(d);

if nargin<2, opts=[]; end
if isfield(opts,'rankeps'),      rankeps      = opts.rankeps;       else rankeps      = numeric_t('1e-12',class_t);  end
if isfield(opts,'showrank'),     showrank     = opts.showrank;      else showrank     =  0;                          end
if isfield(opts,'heuristic'),    heuristic    = opts.heuristic;     else heuristic    =  1;                          end
if isfield(opts,'mingapdif'),    mingapdif    = opts.mingapdif;     else mingapdif    = numeric_t('0.5',class_t);    end
if isfield(opts,'firstlastgap'), firstlastgap = opts.firstlastgap;  else firstlastgap =  0;                          end
if isfield(opts,'fixedrank'),    fixedrank    = opts.fixedrank;     else fixedrank    = -1;                          end
if isfield(opts,'fixedrankgap'), fixedrankgap = opts.fixedrankgap;  else fixedrankgap =  1;                          end

% These options might be set by the calling function for additional info
% that is displayed if showrank > 0 (the names of the calling function and original matrix size)
if isfield(opts,'call'),         call = opts.call;                  else call = '';        end
if isfield(opts,'m1'),           m1 = opts.m1;                      else m1 = 0;           end
if isfield(opts,'m2'),           m2 = opts.m2;                      else m2 = 0;           end

n = length(d);
choice = 0;

% no singular values, rank is 0
if n == 0 
    r = 0; 
    return; 
end;
    
% just one singular value, we compare it to rankeps
if n == 1
    if d > rankeps
        r = 1;
    else
        r = 0;
    end
    return;
end

% zero rank - all singular values are very small, rank is zero
if d(1) < 1000*rankeps 
    r = 0;
    if showrank
        fprintf('(%4d x %4d) r: %4d  gap: %5.1e mgap: %5.1e pos: %4d d(1): %5.1e d(r): %5.1e d(r+1): %5.1e ingap: %5.1e ver: %d avgap: %5.1e %s\n',...
          m1,m2,0,0,0,0,d(1),0,d(1),0,0,0,call);
    end
    return;
end;

h = d/d(1); % normalized singular values

% if option is set to search for the first rank, we lower rankeps so that
% there is no gap between the smallest d(n) and rankeps 
% in such case we can have full rank only if all remaining gaps are small
if heuristic && (firstlastgap == 1) && (h(n)>0) % in this option we don't want full rank 
    if h(n)>rankeps
        rankeps=h(n);
    end
end

% rank is the number of singular values, larger than rankeps*d(1)
% this is the classical choice if heuristic=0 and fixedrank is not used
r = sum(d > rankeps * d(1));
choice = 0;

% we compute the gaps (log10) and set the gaps for singular values below
% rankeps to zero, so that they can not be selected

warning off % turn off log of zero warning in old Matlab versions
gaps = log10(h(1:end-1)./h(2:end));
gaps(n,1) = log10(h(n)/rankeps);
gaps(h<rankeps) = 0;
warning on

[maxgap,posgap] = max(gaps);
avggap  = -log10(rankeps)/n; % average gap

findgap = 1;
if (fixedrank >= 1) && (fixedrank <= n)
    findgap = 0;
    % we select the prescribbed rank unless the gap too small indicates that this is wrong
    if fixedrank == n
        if h(n)<rankeps/100
            findgap = 1;
        end
    else
        if gaps(fixedrank)<fixedrankgap
            findgap = 1;
        end
    end
    if findgap == 0
        r = fixedrank;
        choice = -1;
    else
        if ~heuristic
            % we overrule fixed rank with the classic choice
            choice = 1;
        end
    end
end

if findgap && heuristic
   % if the largest gap is much larger from the second largest, we have a
   % "clear" gap and we take this position for rank
   cleargap = 0;
   choice = 0;
   [tmp,pos] = sort(gaps,'descend');
   if tmp(1)> tmp(2) + mingapdif
       cleargap = 1;
       choice = 2;
       r = pos(1);
   end

   % if gap is not clear, we may use special heuristics if selected in options
   if ~cleargap && firstlastgap
      pos = find(gaps > 1.5*avggap); % positions of gaps that are above average
      if firstlastgap == 1
          % with this option we select the first gap that is above average, if
          % none of them stands out, we say that rank is zero
          if isempty(pos)
             r = 0;
             choice = 3;
          else
             r = pos(1);
             choice = 4;
          end
      end
      if firstlastgap == 2
          % with this option we select the last gap that is above average, if
          % none of them stands out, we say that rank is full
          if isempty(pos)
             r = n;
             choice = 5;
          else
             r = pos(end);
             choice = 6;
          end
      end
   end     
end

% Display diagnostics
if showrank
   spm1 = max(r-2,1);
   zgm1 = min(r+2,length(d));
   interval = unique([1:min(3,length(d)) spm1:zgm1 length(d)]);
   gaps(end+1,1) = numeric_t('Inf',class_t);
   if showrank==2 
      disp([interval' d(interval) d(interval)/d(1) gaps(interval)])
   end
   if (r<length(d)) && (r>0)
       currgap = gaps(r);
   else
       currgap = 0;
   end
   if r==0
       dran = 0;
   else
       dran = d(r);
   end
   ingap = max(gaps(1:end-2));
   if r<n 
       nextd = d(r+1);
   else
       nextd = 0;
   end;
   fprintf('(%4d x %4d) r: %4d  gap: %5.1e mgap: %5.1e pos: %4d d(1): %5.1e d(r): %5.1e d(r+1): %5.1e ingap: %5.1e ver: %d avgap: %5.1e %s\n',...
       m1,m2,r,currgap,maxgap,posgap,d(1),dran,nextd,ingap,choice,avggap,call);
end

end


