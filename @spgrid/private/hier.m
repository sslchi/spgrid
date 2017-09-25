function [ss,ll] = hier(l,type)
% HIER  Return index of grids in Hierarchical order.
%  [SS,LL] = HIER(L,TYPE), where L is the level, SS is the index of grids in
%  Hierarchical order. LL(i) is the level of grids SS(i) on. Note the NO. of
%  total grids are 2^L+1. TYPE is the type of sparse grids, it should be one of
%  {'disall','disinner','quad','symolyak'}.
% 
% Example
%  l = 3, type = 'disinner'
% level
%   0       o                       o
%   1                   o
%   2             o           o
%   3          o     o     o     o
%           o  o  o  o  o  o  o  o  o
%           1  2  3  4  5  6  7  8  9
%
%
% [ss,ll] = spgrid.hier(3,'disinner')


% Checked: 07-Sep-2017.
% $Last revised: 07-Sep-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This algorithm is faster than the one given by Jie Shen and Haijun Yu 
% list below this program. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    
% initialisation
N = 2^l; tt = 2:N; tt = tt(:); 
ss = zeros(N+1,1); ll = ss;

% level zero
ss([1,2]) = [1,N+1];

for k = (l-1):-1:0 
    Nb = 2^(k) + 2; Ne = 2^(k+1) +1;
    ss(Nb:Ne) = tt(1:2:end); tt(1:2:end) = [];   
    ll(Nb:Ne) = (k+1)*ones(Ne-Nb+1,1);  
end


switch type
    case {'disall','disinner'}
        return
    case {'quad','symolyak'}
        if l == 0
            ss = 1; ll = 0;
        else
            ss([1,2,3]) = ss([3,1,2]); ll([1,2,3]) = [0,1,1];
        end
end





end

%% ----------------------------- algorithm of Shen Jie -------------------------
%%%  Ref: Efficient	Spectral Sparse	Grid Methods and Applications to
%%%  High-Dimensional Elliptic Problems. Shen Jie & HaiJun Yu. SIAM J. Sci.
%%%  Comput.

% function s = hier(l)
% 
% N = 2^l;
% 
% s = zeros(N+1,1); s(1) = 1;
% 
% for k = 1:N
%     m = ceil(log2(k));
%     s(k+1) = 2^l/2^m*(1+2*(2^m-k))+1;
% end
%     
% 
% end
