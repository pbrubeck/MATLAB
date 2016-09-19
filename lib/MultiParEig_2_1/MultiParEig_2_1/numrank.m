function r = numrank(d, opts)

%NUMRANK   Find a gap in singular values
%
% r = NUMRANK(d, opts) returns rank from a set of singular values d
%
% Arguments: 
%   d : a set of singular values (e.g., result of the function svd)
%   opts: options
% 
% Options in opts:
%   - rankeps : 1e-12 (discard all singular values smaller from rankeps)
%   - heuristic : 1 (take the largest gap if it is clear) 
%   - mingap : 0.5 (minimum difference between the largest and the second largest gap, where gap(i) = log10(d(i)/d(i+1))
%   - firstlastgap : 0 (special strategy for selecting rank - default is 0)
%       1 : if gap is not clear, select first gap - use when it is safer to detect smaller rank and you know that the rank can not be full
%       2 : if gap is not clear, select last gap - use when it is safer to detect larger rank
%   - showrank : 0 (display diagnostics, 1: partial, 2: detail)
%
% Output:
%   r : numerical rank

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<2, opts=[]; end
if isfield(opts,'rankeps'),      rankeps = opts.rankeps;            else rankeps = 1e-12;  end
if isfield(opts,'showrank'),     showrank = opts.showrank;          else showrank = 0;     end
if isfield(opts,'heuristic'),    heuristic = opts.heuristic;        else heuristic = 1;    end
if isfield(opts,'mingap'),       mingap = opts.mingap;              else mingap = 0.5;     end
if isfield(opts,'firstlastgap'), firstlastgap = opts.firstlastgap;  else firstlastgap = 0; end

% These options migth be set by the calling function for additional info
% that is displayed if showrank > 0 (the names of the calling function and original matrix size)
if isfield(opts,'call'),         call = opts.call;                 else call = '';       end
if isfield(opts,'m1'),           m1 = opts.m1;                     else m1 = 0;          end
if isfield(opts,'m2'),           m2 = opts.m2;                     else m2 = 0;          end

n = length(d);

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
        fprintf('(%4d x %4d) r:%4d  d(1): %5.1e %s\n',m1,m2,0,double(d(1)),call);
    end
    return;
end;

h = d/d(1); % normalized singular values

% if option is set to search for the first rank, we lower rankeps so that
% there is no gap between the smallest d(n) and rankeps 
% in such case we can have full rank only if all remaining gaps are small
if (firstlastgap == 1) && (h(n)>0) % in this option we don't want full rank 
    if h(n)>rankeps
        rankeps=h(n);
    end
end

% rank is the number of singular values, larger than rankeps*d(1)
% this the classical choice if no other options are set
r = sum(d > rankeps * d(1));
choice = 0;

% we compute the gaps (log10) and set the gaps for singular values below
% rankeps to zero, so that they can not be selected

warning off % turn off log of zero warning in old Matlab versions
gaps = log10(h(1:end-1)./h(2:end));
gaps(n,1) = log10(h(n)/rankeps);
gaps(h<rankeps) = 0;
warning on

[maxgap,posgap] = max(double(gaps));
avggap  = -log10(rankeps)/n; % average gap

if heuristic
   % if the largest gap is much larger from the second largest, we have a
   % "clear" gap and we take this position for rank
   cleargap = 0;
   choice = 0;
   [tmp,pos] = sort(gaps,'descend');
   if tmp(1)> tmp(2) + mingap
       cleargap = 1;
       choice = 1;
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
             choice = 2;
          else
             r = pos(1);
             choice = 3;
          end
      end
      if firstlastgap == 2
          % with this option we select the last gap that is above average, if
          % none of them stands out, we say that rank is full
          if isempty(pos)
             r = n;
             choice = 4;
          else
             r = pos(end);
             choice = 5;
          end
      end
   end     
end

% Display diagnostics
if showrank
   spm1 = max(r-2,1);
   zgm1 = min(r+2,length(d));
   interval = unique([1:min(3,length(d)) spm1:zgm1 length(d)]);
   gaps(end+1,1) = Inf;
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
       dran = double(d(r));
   end
   ingap = max(double(gaps(1:end-2)));
   if r<n 
       nextd = double(d(r+1));
   else
       nextd = 0;
   end;
   fprintf('(%4d x %4d) r: %4d  gap: %5.1e mgap: %5.1e pos: %3d d(1): %5.1e d(r): %5.1e d(r+1): %5.1e ingap: %5.1e ver: %d avgap: %5.1e %s\n',...
       m1,m2,r,double(currgap),double(maxgap),posgap,double(d(1)),dran,nextd,double(ingap), choice, double(avggap),call)
end



