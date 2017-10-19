function w = findtext(string)
%FINDTEXT   find out files containing specific string
%
% W = FINDTEXT(STRING) will print all the files cointaining
% string STRING, the input must be a string.
% 

% Written by G. Wang, Dec. 2015


if ischar(string)
    a = ['grep -rn ''',string,''' ./ ' ];
    [~,w] = system(a);
else
    error('input should be a string')
end

end