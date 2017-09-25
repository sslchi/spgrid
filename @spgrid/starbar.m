function indset = starbar(d, l, type)
%STARBAR Returns the multi-index set of level L
%   INDSET = STARBAR(D, L, TYPE)
%   This program use star and bars method for levelset of level = L. 

% Checked: 12-Sep-2017
% $Last revised: 12-Sep-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com


% Parse inputs
type  = lower(type);
valid_type = {'disall','disinner','quad','smolyak'};

if ( ~ismember(type,valid_type) )
    error('spgrid:starbar:InPut',...
        'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
end


k = d - 1; n = d + l - 1; v = 1:n;
if k == 0
    indset = l;
else
    indset = nchoosek(v, k);
    m = size(indset, 1);
    indset = [indset,(n+1)*ones(m,1)];
    indset(:,2:end) = indset(:,2:end) - indset(:,1:end-1);
end

% special for each type
if strcmpi(type,'disall')
    indset = spgrid.inner2all(indset);
end

indset = sortrows(indset);
end