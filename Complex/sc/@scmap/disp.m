function disp(f)

% $Id: disp.m 298 2009-09-15 14:36:37Z driscoll $ 

s = char(f);
if isstr(s)
  disp(s)
elseif iscell(s)
  fprintf('\n  SC %s:\n\n',class(f));
  for n = 1:length(s)
    disp(s{n})
  end
end
