function object = numeric_t(expression, class_t)
%NUMERIC_T Creates object of type 'class_t' from 'expression'.
% 
% Input:
%   - expression: numerical entity (scalar, matrix, array) or mathematical 
%                 expression written as string.
%   - classname: target type name, e.g. 'double', 'single', 'mp', etc.
%
% Output: 
%   - object: entity converted to target classname type.
%
% The function allows programmer to write code independent from particular
% numeric type. If code relies on 'numeric_t' then it can be run without changes with 'double', 'single' 
% or even with extended precision type 'mp'. Only 'class_t' parameter
% should be changed, not the code itself.
% 
% Examples:
%
% >> numeric_t('sqrt(2)/2','double')
% ans =
%         0.707106781186548
%
% >> numeric_t('sqrt(2)/2','mp')
% ans = 
%    0.707106781186547524400844362104849
% 
% See other functions in toolbox for more detailed examples.
%
% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

    narginchk(2,2);
    if (strcmpi(class_t,'mp')), object = mp(expression);
    else
        % double & single
        if isnumeric(expression)
            object = expression;
        else
            object = eval(expression);
        end;
    end;
end % numeric_t

