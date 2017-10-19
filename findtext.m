function w = findtext(string)
%FINDTEXT   Find out files containing specific string
%
%   W = FINDTEXT(STRING) will print all the files cointaining
%   string STRING, the input must be a char.
% 
% Examples:
%   findtext('findtext(string)')
%   findtext findtext

% Written by G. Wang, Dec. 2015


if ischar(string)
    a = ['grep -rn ''',string,''' ./ ' ];
    [~,w] = system(a);
else
    error('InPut should be a char.')
end

end