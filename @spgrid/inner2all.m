function indset = inner2all(indset)
% INNER2ALL Turn indset for disinner to indset for disall.
%    INDSET = INNER2ALL(INDSET)

[len,d] = size(indset);

indicator = (indset == 1);

more = sum(2.^sum(indicator,2) - 1);

indset_all = zeros(len+more,d);
indset_all(1:len,:) = indset;

for ii = 1:d
    ff = (indset_all(:,ii) == 1);
    Nii = sum(ff); 
    indii = indset_all(ff,:); indii(:,ii) = 0;
    indset_all(len+1:len+Nii,:) = indii;
    len = len + Nii;
end

indset = indset_all;

end