function status = is_numeric_type_supported(classname)
%IS_NUMERIC_TYPE_SUPPORTED Checks if required numeric type is supported by MultiParEig
%
% Input:
%   - classname: type name to check support for, e.g. 'double', 'single', 'mp', etc.
%
% Output: 
%   - status: true if supported, throws exception otherwise.
%
% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

    status = false;
    if strcmpi(classname,'mp') % Check if MCT is installed.
        installed_toolboxes = ver;
        search_results = arrayfun(@(x)all(strcmp(x.Name,'Advanpix Multiprecision Computing Toolbox')),installed_toolboxes);
        status = any(search_results);
        if ~status
            error('\n%s\n%s\n','Numeric type ''mp'' requires Advanpix Multiprecision Computing Toolbox.',...
                               'Please install the toolbox, add it to the path and try again.');
        end;
    elseif strcmpi(classname,'double') || strcmpi(classname,'single')
        status = true;
    end;
end

