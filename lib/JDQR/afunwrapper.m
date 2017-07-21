function [out1]=afunwrapper(x, flag, varargin)
global afun n pfun;
if nargin <2
    out0=cellfun(afun, num2cell(x, 1), 'UniformOutput', false);
    out1=cat(ndims(x), out0{:});
else
    switch(flag)
        case 'set'
            afun=varargin{1};
            n=varargin{2};
            pfun=varargin{3};
        case 'dimension'
            out1=n;
        case 'preconditioner'
            out0=cellfun(pfun, num2cell(x, 1), 'UniformOutput', false);
            out1=cat(ndims(x), out0{:});
        case 'L'
            out1=[];
        case 'U'
            out1=[];
        otherwise
            error(['Unknown flag ''' flag '''.']);
   end
end
end 

