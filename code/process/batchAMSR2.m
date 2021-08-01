%% Batch download AMSR2 sea ice concentrations from the Univeristy of Bremen
%Modified from code from S.D.Brenner, 2020

%% Clean workspace
clearvars -except rootPath profiles wvdata metData
close all;
saveDir = [rootPath, 'data/AMSR2/'];
if ~exist(saveDir, 'dir'); mkdir(saveDir); end 
months = {'jul', 'aug', 'sep', 'oct', 'nov'}; %specify months from which to download data

%% Open FTP link
% (note: if connect through a firewall, there may be issues)
ftpLink = 'seaice.uni-bremen.de';
ftpName = 'anonymous';
ftpPass = '';
ftpObj = ftp(ftpLink, ftpName, ftpPass, 'LocalDataConnectionMethod', 'active');

%% Download latitude/longitude grid
cd(ftpObj, '~/grid_coordinates/n3125/')
mget(ftpObj,'LongitudeLatitudeGrid-n3125-ChukchiBeaufort.hdf', saveDir);

%% Get 2018 files
dataDir2018 = '~/amsr2/asi_daygrid_swath/n3125/2018/'; %Download 2018 data
locDir = '/ChukchiBeaufort/';
saveDir = [saveDir, 'AMSR2raw/'];
if ~exist(saveDir, 'dir'); mkdir(saveDir); end 

%%
for n = 1:length(months)
    dataDir = [dataDir2018, months{n}, locDir];
    cd(ftpObj, dataDir);
    
    files = dir(ftpObj,'*-v5.hdf');
    
    fprintf('Downloading %g files for %s, 2018\n', length(files), months{n});   
    mget(ftpObj,'*-v5.hdf', saveDir);
end
