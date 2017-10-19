function numsp = numsp(levelset, type, level)
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
    
    switch type
        case {'quad','smolayak'}
            level_1 = bsxfun(@eq, levelset, 1);
            level_2 = bsxfun(@eq, levelset, 2);
            levelset(level_1) = 2;
            levelset(level_2) = 3;
            numsp = bsxfun(@power, 2, bsxfun(@minus,levelset, 2) );
        case 'disall'
            level_0 = bsxfun(@eq, levelset, 0);
            levelset(level_0) = 2;
            numsp = bsxfun(@power, 2, bsxfun(@minus,levelset, 1) );
        case 'disinner'
            numsp = bsxfun(@power, 2, bsxfun(@minus,levelset, 1) );           
    end
    numsp_ii = prod(numsp,2);
    numsp = numsp_ii;
    for ii = 2:length(numsp)
        numsp(ii) = numsp_ii(ii) + numsp(ii-1);
    end
    numsp = [0;numsp];
end



if nargin == 3  
    l = type; d = levelset; type = lower(level);
    
    valid_type = {'disinner','disall','quad','smolyak'};
    if ( ~ismember(type, valid_type))
        error('spgrid:numsp:InPut',...
            'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
    end
    
    levelset = spgrid.levelset(d,l,type);
    numsp = spgrid.numsp(levelset,type); numsp = numsp(end);
end


end