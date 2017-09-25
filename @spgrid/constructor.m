function II = constructor(d, l, type)
% CONSTRUCTOR  The main constructor of SPGRID.
%   II = CONSTRUCTOR(D, L, TYPE) returns the index of sparse grids in
%   Hierarchical order.

% Checked: 23-Sep-2017
% $Last revised: 23-Sep-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com

levelset = spgrid.levelset(d, l, type); % set of levels
numsp = spgrid.numsp(levelset,type);    % cumulative number of points for indset.
II = zeros(numsp(end),d);               % initialisation

for ii = 1:size(levelset,1)
    indx = levelset(ii,:);
    yy = cell(d,1);
    for kk = 1:d
        ll = indx(kk); 
        switch type
            case {'quad','symolyak'}
                if (ll == 1)
                    yy{kk} = 1;
                elseif (ll == 2)
                    yy{kk} = [2;3];
                else
                    yytemp = 2^(ll-2):(2^(ll-1) - 1);
                    yy{kk} = bsxfun(@plus, yytemp(:), 2);
                end
            case {'disall','disinner'}
                if (ll == 0)
                    yy{kk} = [1;2];
                elseif (ll == 1)
                    yy{kk} = 3;
                else 
                    yytemp = 2^(ll-1):(2^ll - 1);
                    yy{kk} = bsxfun(@plus, yytemp(:), 2);
                end
                
        end
                
    end
    nodes = yy{1} ;
    for kk = 2:d
        newnodes = yy{kk};
        nodes = [kron(nodes,ones(size(newnodes,1),1)),...
            kron(ones(size(nodes,1),1),newnodes)];
    end    
    temp_a = numsp(ii) + 1; temp_b = numsp(ii+1);
    II(temp_a:temp_b,:) = nodes;
end


end