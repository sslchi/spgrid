function [type, number] = parsespin(spvar, par)
% PARSESPIN Parse input of SPGRID functions.
%   TYPE = PARSESPIN( SPVAR, PAR) will return the result of varargin.
%   SPVAR should be a cell structure and all the arry of SPVAR should be char or
%   string, PAR is a char or string.
%
% NOTE:
%   The defalut values are also set by this function.

%% Parse input

if ischar(par)
    par = lower(par);
else
    error('SPGRID:parsespin:InPut',...
        'Second input should be CHAR.')
end

if ( ~strcmpi(par,'allchar') )
    if iscellstr(spvar)
        spvar = lower(spvar);
    else
        for i = 1:length(spvar)
            if (~ischar(spvar{i}))
                if ( isnumeric(spvar{i}) )
                    spvar{i} = num2str( spvar{i} );
                else
                    spvar(i) = [];
                    warning('qmesh:square',...
                        ['Inputs after the second input',...
                        'should be a CHAR or NUMERIC, otherwise, it will be ignored.']);
                end
            end
        end
        spvar = lower(spvar);
    end
end


%% switch case
switch par
    %% Determine the type of sparse grid
    case 'type'
        if iscellstr(spvar)
            spvar = lower(spvar);
        else
            error('spgrid:input',...
                'Input(varargin) for spgrid contains variable not a char.\n')
        end
        type = 'disall';
        valid_type = {'disall', 'quad','disinner','smolyak'};% get the type of points
        itype = ismember(valid_type, spvar);
        
        if (sum(itype) == 1)
            type = valid_type{itype}; % get domain
        elseif (sum(itype) > 1)
            error('spgrid:constructor',...
                'Please give only one type.')
        end
        if strcmpi(type,'smolyak')
            type = 'quad';
        end
        number = [];
        %% Determine whether all arrays of cell are char and turn numeric to char.
    case 'allchar'
        if iscellstr(spvar)
            type = lower(spvar); number = true;
        else
            isallchar = true;
            for i = 1:length(spvar)
                if (~ischar(spvar{i}))
                    if ( isnumeric(spvar{i}) )
                        spvar{i} = num2str( spvar{i} );
                    else
                        spvar{i} = []; isallchar = false;
                        warning('spgrid:input',...
                            ['Input(varargin) for SPGRID contains variables that ',...
                            'are neither CHAR nor NUMERIC.\n'])
                    end
                end
            end
            type = lower(spvar); number = isallchar;
        end
        
    otherwise
        error('spgrid:parsespin','second input must be a CHAR type/allchar...\n')
end

end