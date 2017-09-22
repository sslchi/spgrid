function indset = indx( d, l, type )
%INDX Give level-index of  sparse grid.
%
% Checked: 12-Sep-2017
% $Last revised: 12-Sep-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com


%  disall & disinner:  D <= |level| <= L + D - 1.
%                --->  indset := level
%                --->  |indset  - 1| < = L - 1   
%                --->  |indset|_inf = L    |indset|_min = 1
%
%  quad & Smolyak :    D <= |level| <= L + D 
%                --->  |level-1| < = L 
%                --->  indset := level-1
%                --->  |indset|_inf = L   |indset|_min = 0;
%
%
% for Smolyak, the disadjoint points in each level is like this
% level                                     l = level - 1.
%   2       o                       o       1   include
%   1                   o                   0
%   3             o           o             2
%   4          o     o     o     o          3
%           o  o  o  o  o  o  o  o  o
%           1  2  3  4  5  6  7  8  9
%
% for disinner, the disadjoint points in each level is like this
% level                                       l = level
%   0       o                       o         0   exclude
%   1                   o                     1
%   2             o           o               2
%   3          o     o     o     o            3
%           o  o  o  o  o  o  o  o  o
%           1  2  3  4  5  6  7  8  9

% Parse inputs
type  = lower(type);
valid_type = {'disall','disinner','quad','smolyak'};

if ( ~ismember(type,valid_type) )
    error('spgrid:indx:InPut',...
        'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
end

indset = [];

if ismember(type,{'disinner','disall'})
    for s = 0:(l-1)
        index = spgrid.starbar(d,s,type);
        indset = [indset;index];
    end
elseif ismember(type,{'quad','smolyak'})
    for s = 0:l
        index = spgrid.starbar(d,s,type);
        indset = [indset;index];
    end
end


end
