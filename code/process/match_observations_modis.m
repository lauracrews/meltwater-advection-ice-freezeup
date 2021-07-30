%Compares in situ sea surface temperature to MODIS SST averaged in a box
%surrounding the observation. Iterates through daily MODIS images to find 
%the closest MODIS SST in time to each observation, up to a threshold number 
%of days away from the observation day. 

close all; clc
clearvars -except rootPath AMSR2 profiles wvdata metData
plotPlatformComparison = false;

%Empty data arrays to hold all of the observations
lons = []; lats = []; times = []; CTs = []; temps = []; platform = [];
startTime = datenum('sept 19 2018')+.3; 
endTime = datenum('oct 15 2018')-.4;
    
timeThreshold = 3; %Number of days within which the MODIS observation must be to the observation
dayOffset = [0:1:timeThreshold; 0:-1:-timeThreshold];
dayOffset = dayOffset(:); dayOffset = dayOffset(2:end); %Array of [0 1 -1 2 -2 ... n -n] where n is timeThreshold. Will iterate through later

%Load underway data from the Healy system 
% metData = load_correct_underwayData;
healyIncrement = 15; %Only use every 15th Healy measurement, so that the number of measurements is comparable across platforms

%Restrict Healy  data to that collected in the time limits
[~, healyStartInd] = min(abs(metData.times - startTime));
[~, healyEndInd] = min(abs(metData.times - endTime));

%Add Healy data to the growing arrays holding data from all platforms
lons = [lons; metData.lons(healyStartInd:healyIncrement:healyEndInd)]; lats = [lats; metData.lats(healyStartInd:healyIncrement:healyEndInd)];
times = [times; metData.times(healyStartInd:healyIncrement:healyEndInd)];
CTs = [CTs; metData.CTs(healyStartInd:healyIncrement:healyEndInd)]; temps = [temps; metData.temps(healyStartInd:healyIncrement:healyEndInd)];
platform = [platform; ones(length(healyStartInd:healyIncrement:healyEndInd), 1)];

%Seaglider and uCTD data
% profiles = loadProfiles;

%Restrict profiles to those within the time range
inTimeMask = zeros(size(profiles.times));
inTimeMask(profiles.times >= startTime & profiles.times < endTime) = 1;
profNums = find(inTimeMask == 1 & profiles.qualFlag == 1);

%Add Seaglider and uCTD data to the arrays holding data from all platforms
lons = [lons; profiles.lons(profNums)]; lats = [lats; profiles.lats(profNums)]; 
times = [times; profiles.times(profNums)];
CTs = [CTs; nanmean(profiles.CT(1:3, profNums))']; temps = [temps; nanmean(profiles.temps(1:3, profNums))'];
platform = [platform; 2 .* ones(size(profiles.lons(profNums)))];

%Wave Glider data
% wvdata = loadWaveglider;
[SA, ~, ~] = gsw_SA_Sstar_from_SP(wvdata.salts(:, 3), gsw_p_from_z(-wvdata.depths(:,  3), wvdata.lats), wvdata.lons, wvdata.lats);
CT = gsw_CT_from_t(SA, wvdata.temps(:,  3), gsw_p_from_z(-wvdata.depths(:,  3), wvdata.lats));

%Add Wave Glider data to the arrays holding data from all platforms
lons = [lons; wvdata.lons]; lats = [lats; wvdata.lats];
times = [times; wvdata.times];
temps = [temps; wvdata.temps(:, 1)]; %At surface
CTs = [CTs; CT]; %at 9 m (salinity required to calculate CT) so not used
platform = [platform; 3 .* ones(size(wvdata.lons))];

%Load daily SST files
files_daily = dir([rootPath, 'data/ModisTerra/daily/*DAY*.nc']);
files_daily  = {files_daily.name};
if isempty(files_daily)
     disp(['MODIS-Terra sea surface temperature data not found. Download daily files for 19 Sept 2018 to 16 Oct 2018 from', ...
            newline, 'https://podaac-tools.jpl.nasa.gov/drive/files/allData/modis/L3/terra/11um/v2019.0/4km/daily/2018', ...
            newline, 'to the directory ~/meltwaterAdvection/data/ModisTerra/daily/'])
        
    return
end

%Iterate through files and extract the date for each by parsing the file name
modisDays = nan .* ones(size(files_daily)); %Convert time of each image
for i = 1:length(files_daily)
    curFile  =  files_daily{i};
    modisDays(i) = datenum([curFile(17:18), ' ', curFile(19:20), ' ', curFile(13:16)]);
end

%% Iterate through all observations, find closest MODIS SST to each observation (if a MODIS observation exists)
clim_temp = [-1.5, 1]; %color limits for plotting temperatures

modisSST_nearObservations = nan .* ones(size(lons));
modisTime_nearObservations = nan .* ones(size(lons));
makePlot = false; close all
for i = 1:length(lons)
    disp(['Matching observation ', num2str(i), ' of ', num2str(length(lons)), ' to MODIS SST'])
    
    %This is the geographic range around each point to search for 
    %MODIS data to average. Box is roughly 15 km
    minlon = lons(i) - 0.25; maxlon = lons(i) + 0.25;  %At 75 N, 0.15 degrees latitude ~= 16 km, 0.5 degrees longitude ~= 14.5 km
    minlat = lats(i) - 0.075; maxlat = lats(i) + 0.075; 
    pt1 = [minlon, minlat]; pt2 = [maxlon, minlat]; pt3 = [maxlon, maxlat]; pt4 = [minlon, maxlat];
    pts = [pt1(1), pt2(1), pt3(1), pt4(1); pt1(2), pt2(2), pt3(2), pt4(2)];
       
    if makePlot; m_proj('lambert', 'lon', [minlon-1 maxlon+1], 'lat', [minlat-.25 maxlat+.25]); end
   
    %Search from closest MODIS day in time to more distant in time MODIS days 
    %On each day, load and average MODIS data within pts - if there is data
    %in that range (so the average is not nan) then leave the loop. 
    modisMean = nan;
    count = 1;
    while count <= length(dayOffset) & isnan(modisMean) 
        
        %Identify which MODIS file to load, load the file, identify points
        %within the geographic bounds
        [~, modisInd] = min(abs(modisDays - (times(i)+dayOffset(count))));
        curFile = files_daily{modisInd};
        [modisLon, modisLat, modisSST, ~, ~, modisTitle] = loadMODISsst(curFile, minlon-.5, maxlon+.5, minlat-.125, maxlat+.125);
        [in, on] = inpolygon(modisLon, modisLat, pts(1, :), pts(2, :));
        inBoxMask = in;

        %Average MODIS data within the bounds  - if this is not nan (so
        %there is some data) then we are done and will break out of the loop
        modisMean = nanmean(modisSST(inBoxMask));
        
        %Optionally make plots showing what MODIS data is available in the vicinity
        if makePlot
            subplot(2, 4, count)
            h1 = m_pcolor(modisLon, modisLat, modisSST);
            hold on  
            m_contour(modisLon, modisLat, double(inBoxMask), 'k')
            caxis(clim_temp); cmocean('thermal')
            freezeColors(h1)
            m_grid('linestyle', 'none', 'fontsize', 12) %Need to m_grid on both sets of axes 
            m_scatter(lons(i), lats(i), 60, 'w', 'filled')
            m_scatter(lons(i), lats(i), 40, CTs(i), 'filled')
            title([modisTitle, sprintf('\n'), 'Nearby SST = ', num2str(modisMean)])            
        end
        count = count + 1;
    end
    
    %Record the day of the MODIS data matched to the observation as well as
    %the average MODIS data near the observation
    modisSST_nearObservations(i) = modisMean;
    modisTime_nearObservations(i) = modisDays(modisInd);
        
end

%This takes a little while for all of the observations - save to reload later for speed
cd([rootPath, 'data/'])
save('modisComparison.mat', 'lons', 'lats', 'times', 'temps', 'CTs', 'platform', 'startTime', 'endTime', 'timeThreshold', 'modisSST_nearObservations', 'modisTime_nearObservations');
