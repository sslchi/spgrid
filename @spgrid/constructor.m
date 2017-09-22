function II = constructor( d, l, type )
% CONSTRUCTOR  The main constructor of SPGRID.
%   II = CONSTRUCTOR( D, L,TYPE)
%   The multi_index, i.e., level, satisfy D <= |level|_1 <= L + D 

% Checked: 12-Sep-2017
% $Last revised: 12-Sep-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com


[ss,ll] = spgrid.hier(l,type);
indset = spgrid.indx(d,l,type);% for set of multi_index
numsp = spgrid.numsp(indset,type); % cumulative number of points for indset.
II = zeros(numsp(end),d);

for ii = 1:size(indset,1)
    indx = indset(ii,:);
    indicator = bsxfun(@eq, ll, indx);
    
    yy = cell(d,1);
    for kk = 1:d
        yy{kk} = ss(indicator(:,kk));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    nodes = yy{1} ;
    for kk = 2:d
        newnodes = yy{kk};
        nodes = [kron(nodes,ones(size(newnodes,1),1)),...
            kron(ones(size(nodes,1),1),newnodes)];
    end 
%%%-------- equals to the following, but the above is more efficient. ----------
%     
%    
%         decim = sum(indicator,1);
%     
%         Nii = prod(decim);
%         nodes = zeros(Nii,d);
%     
%        for kk = 1:Nii
%            for tt = d:-1:1
%     
%               ptt = decim(tt);
%     
%               if (tt == d)
%                   cptt = 1;
%               else
%                   cptt = prod(decim((tt+1:end)));
%               end
%     
%               idii = mod(floor((kk-1)/cptt),ptt) + 1;
%               yytemp = yy{tt};
%               nodes(kk,tt) = yytemp(idii);
%            end
%        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    temp_a = numsp(ii) + 1; temp_b = numsp(ii+1);
    II(temp_a:temp_b,:) = nodes;
end


end