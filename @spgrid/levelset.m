function levelset = levelset( d, l, type )
%LEVELSET Give set of levels for  sparse grid.
%
% Checked: 23-Sep-2017
% $Last revised: 23-Sep-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com


% Parse inputs
type  = lower(type);
valid_type = {'disall','disinner','quad','smolyak'};

if ( ~ismember(type,valid_type) )
    error('spgrid:indx:InPut',...
        'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
end

levelset = [];


for s = 0:l
    index = spgrid.starbar(d,s,type);
    levelset = [levelset;index];
end

end
