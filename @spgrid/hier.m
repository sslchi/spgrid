function [hid,level] = hier(l,type)
% HIER  Return index of grids in Hierarchical order.
%  [HID,LEVEL] = HIER(L,TYPE), where L is the level, HID is the index of grids in
%  Hierarchical order. LEVEL(i) is the level of grids HID(i) on. TYPE is the
%  type of sparse grids, it should be one of
%  {'disall','disinner','quad','symolyak'}.
%  X(hid) will arrange Nodes in Hierarchical order, where X is the nodes in
%  ascending order.



% Checked: 24-Sep-2017.
% $Last revised: 24-Sep-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com

%% For disinner & disall:
%                      D <= |levelset|_1 <= L + D
%                --->  |levelset  - 1|_1 < = L    }
%                      |levelset|_min = 1         }
%                --->  |levelset|_inf = L + 1
%
% the disadjoint points in each level is like this:
%
%                                           _____levelset_____
%   L                                       disinner    disall
%           o                       o                      0
%   0                   o                       1          1
%   1             o           o                 2          2
%   2          o     o     o     o              3          3
%           o  o  o  o  o  o  o  o  o
%           1  6  4  7  3  8  5  9  2       Nodes in Hierarchical order
%           1  2  3  4  5  6  7  8  9       Nodes in oordinary order
%
%
%% For quad & Smolyak :
%                      D <= |levelset|_1 <= L + D
%                --->  |levelset - 1|_1 < = L   }
%                      |levelset - 1|_min = 0   }
%                --->  |levelset - 1|_inf = L
%
% the disadjoint points in each level is like this:
%
%   L                                       levelset
%
%   1       o                       o          2
%   0                   o                      1
%   2             o           o                3
%   3          o     o     o     o             4
%           o  o  o  o  o  o  o  o  o
%           2  6  4  7  1  8  5  9  3       Nodes in Hierarchical order
%           1  2  3  4  5  6  7  8  9       Nodes in oordinary order

% Parse input
if l < 0
    error('spgrid:hier:input','L must be a no-negative integer.')
end

type = lower(type);
valid_type = {'disall','disinner','quad','symolyak'};
if ( ~ismember(type, valid_type) )
    error('spgrid:hier:InPut',...
        'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
end
    


if ismember(type, {'disall','disinner'})
    [hid, level] = sub_hier(l+1);
elseif ismember(type, {'quad','symolyak'})
    [hid, level] = sub_hier(l);
    level = bsxfun(@plus, level, 1);
    if (l == 0)
        hid = 1; level = 1;
    else
        hid([1,2,3]) = hid([3,1,2]);
        level([1,2,3]) = [1,2,2];
    end
end

end


function [hid, level] = sub_hier(l)
%SUB_HIER  return hid  Hierarchical order and their level.

N = 2^l;
tt = 2:N; tt = tt(:);
hid = zeros(N+1,1);
level = zeros(N+1,1);
hid([1,2]) = [1,N+1];
for k = (l-1):-1:0 % loop from the last level
    Nb = 2^(k) + 2; Ne = 2^(k+1) +1;
    hid(Nb:Ne) = tt(1:2:end); tt(1:2:end) = [];
    level(Nb:Ne) = (k+1)*ones(Ne-Nb+1,1);
end

end

