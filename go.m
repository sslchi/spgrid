function go(dir)
%GO Go to specified directory.
%   GO('home') go to the home  directory.
%   GO(dir) go to the directory dir, where dir subfolder of ESFEM.
%
% Examples:
%   GO data 
%   or use
%   GO('data')


% $Date: 2017/03/24$
% Copyright (c) G. Wang, Email: wangguanjie0@126.com


cd('/home/dell/Documents/github/distrib')

if nargin == 0
    go home; return;
end

if strcmpi(dir,'home'); return; end

if ~isdir(['./',dir])
    go home; return;
end

cd(dir)

end
