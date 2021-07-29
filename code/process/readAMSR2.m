%% Read AMSR2 sea ice concentraions
% S.D.Brenner, 2020. Minor adaptations by Laura Crews

%% Clean workspace;
clearvars; clc;
makePlots = false;
saveName = 'AMSR2_2018.mat';
saveDir = [rootPath, 'data/'];

%% Load latitude/longitude grid data
AMSR2dir = [rootPath, 'data/AMSR2/'];
cd(AMSR2dir)

gridFName = 'LongitudeLatitudeGrid-n3125-ChukchiBeaufort.hdf';
lat = hdfread(gridFName,'Latitudes');
lon = hdfread(gridFName,'Longitudes');

%% Loop and read in sea ice concentration data

dataDir = [rootPath, 'data/AMSR2/AMSR2raw/'];
files = dir([dataDir,'*.hdf']);
numFiles = length(files);

% Pre-allocate
mattime = NaN(1,numFiles);
SIC = NaN([size(lat),numFiles]);

% Loop and read
for n = 1:numFiles
	fName = files(n).name;
    nSIC = hdfread([dataDir,fName],'ASI Ice Concentration')/100;
    SIC(:,:,n) = nSIC;
    
    mattime(n) = datenum( fName(17:24) , 'yyyymmdd' );
end

%% Make and save structure
AMSR2.lat = lat;
AMSR2.lon = lon;
AMSR2.mattime = mattime;
AMSR2.SIC = SIC;

cd(saveDir)
save(saveName,'AMSR2');

%% Visualize
if makePlots
    close all
    m_proj('lambert','lon',[-187.5,-122.5],'lat',[69,77.5],'rect','on');   
    CL = [ -0.0100, 1 ];
    for n = 1:numFiles
        clf
        nSIC = SIC(:,:,n);
        m_pcolor(lon, lat, nSIC, 'linestyle', 'none');
        set(gca,'clim',CL);
        cmocean('ice')
        m_grid;
        cb = colorbar;
        ylabel(cb,'Sea ice concentration');
        
        drawnow();
        pause(0.1)
    end
end

