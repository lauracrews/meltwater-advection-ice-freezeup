%pwp_ReadOutput Routine to read output of Matlab version of PWP
%
%
% Open Binary Files
fidu = fopen([runname '_u'],'r','b');
fidv = fopen([runname '_v'],'r','b');
fidt = fopen([runname '_t'],'r','b');
fids = fopen([runname '_s'],'r','b');

% Load Data
V = fread(fidv,[nz,inf],'float32');
U = fread(fidu,[nz,inf],'float32');
temp = fread(fidt,[nz,inf],'float32');
sal  = fread(fids,[nz,inf],'float32');

% Close them all
fclose('all');
clear fid*

delete([runname '_u'])
delete([runname '_v'])
delete([runname '_t'])
delete([runname '_s'])

% Write MAT File
%save([runname '_Output'],'U','V','time','z','temp','sal');
save([runname '_Output'],'U','V','time','z','dz','dt','temp','sal','config','forc', 'init', 'zmld', 'zmld_static', 'runname','nn*');
