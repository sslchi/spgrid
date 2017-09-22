function numsp = numsp( indset, type, level )
% NUMSP give the number of sparse points.
%   NUMSP = NUMSP( INDSET ), here INDSET is the set of multi_index for dimension D
%   level L sparse grids. NUMSP will give the number of sparse grids for each
%   multi_index, i.e., INDSET(i,:).
% 
%   NUMSP = NUMSP( D, L, TYPE) will give the total number of sparse grids for
%   dimension D, level L and TYPE sparse grid. TYPE must be one of
%   {'disinner','disall','quad','smolyak'}.

% Checked: 11-Sep-2017
% $Last revised: 11-Sep-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com



if nargin == 2
    type = lower( type ); 
    valid_type = {'disinner','disall','quad','smolyak'};
    if ( ~ismember(type, valid_type))
        error('spgrid:numsp:InPut',...
            'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
    end
    if strcmpi(type,'disall')
        ind_0 = (indset == 0);
        indset(ind_0) = 2;
        numsp = 2.^(indset - 1);
    elseif strcmpi(type,'disinner')
        numsp = 2.^(indset -1);
        
    elseif ismember(type,{'quad','smolyak'})
        ind_0 = (indset == 0);
        ind_1 = (indset == 1);
        indset(ind_0) = 1;
        indset(ind_1) = 2;
        numsp = 2.^(indset - 1);
        
    end
    numsp_ii = prod(numsp,2);
    numsp = numsp_ii;
    for ii = 2:length(numsp)
        numsp(ii) = numsp_ii(ii) + numsp(ii-1);
    end
    numsp = [0;numsp];
end



if nargin == 3  
    l = type + 1; d = indset; type = lower(level);
    
    valid_type = {'disinner','disall','quad','smolyak'};
    if ( ~ismember(type, valid_type))
        error('spgrid:numsp:InPut',...
            'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
    end
    
    indset = spgrid.indx(d,l,type);
    numsp = spgrid.numsp(indset,type); numsp = numsp(end);
end


end