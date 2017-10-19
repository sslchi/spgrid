function setup
%SETUP Add all m-files of this project to MATLAB search path.

% path and set the home dir in program go.m
cdir = fileparts(which('setup.m'));
% cdir2 = fileparts(which('guide.pdf'));
% 
% if strcmp(cdir1,cdir2)
%     cdir = cdir1;
% else
%     fprintf('Can not get the home directory, set up by hands.\n.');
% end

cd(cdir);



%set the go.m function
fid = fopen('go.m','w');


        

fprintf(fid, 'function go(dir)\n');
fprintf(fid,'%%GO Go to specified directory.\n');
fprintf(fid,'%%   GO(''home'') go to the home  directory.\n');
fprintf(fid,'%%   GO(dir) go to the directory dir, where dir subfolder of ESFEM.\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% Examples:\n');
fprintf(fid,'%%   GO data \n');
fprintf(fid,'%%   or use\n');
fprintf(fid,'%%   GO(''data'')\n\n\n');

fprintf(fid,'%% $Date: 2017/03/24$\n');
fprintf(fid,'%% Copyright (c) G. Wang, Email: wangguanjie0@126.com\n\n\n');


fprintf(fid,'cd(''%s'')\n\n',cdir);

fprintf(fid,'if nargin == 0\n');
fprintf(fid,'    go home; return;\n');
fprintf(fid,'end\n\n');

fprintf(fid,'if strcmpi(dir,''home''); return; end\n\n' );

fprintf(fid,'if ~isdir([''./'',dir])\n');
fprintf(fid,'    go home; return;\n');
fprintf(fid,'end\n\n');



fprintf(fid,'cd(dir)\n\n');
fprintf(fid,'end\n');

fclose(fid);

% set path
cd(cdir)

addpath(genpath(pwd));

savepath
clear 
