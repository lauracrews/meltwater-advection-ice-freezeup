%Makes plot of SMOS sea surface salinity along with sea surface
%salinity from in situ observations. Iterates through weekly SMOS images and finds 
%in situ observations taken at the time of those images. Options to plot data from all 
%platforms on one figure or to make separate panels for Healy underway data, 
%Seaglider and uCTD data, and Wave Glider data. 
clearvars -except rootPath AMSR2 profiles wvdata metData
plotPlatformComparison = false;

saveFigs = true;
saveDir = [rootPath, 'figures/fig3/'];
saveName = 'satellite_inSitu_comparison';

%Empty data arrays to hold all of the observations
lons = []; lats = []; times = []; SAs = []; salts = []; platform = [];
startTime = datenum('sept 19 2018')+.3; 
endTime = datenum('oct 15 2018')-.4;
    
timeThreshold = 0; %Number of days within which the SMOS observation must be to the observation
defineSODAconstants

%Load underway data from the Healy system 
% metData = load_correct_underwayData;
healyIncrement = 15; %Only use every 15th Healy measurement, so that the number of measurements is comparable across platforms

%Restrict Healy  data to that collected in the time limits
[~, healyStartInd] = min(abs(metData.times - startTime));
[~, healyEndInd] = min(abs(metData.times - endTime));

%Add Healy data to the growing arrays holding data from all platforms
lons = [lons; metData.lons(healyStartInd:healyIncrement:healyEndInd)]; lats = [lats; metData.lats(healyStartInd:healyIncrement:healyEndInd)];
times = [times; metData.times(healyStartInd:healyIncrement:healyEndInd)];
SAs = [SAs; metData.SAs(healyStartInd:healyIncrement:healyEndInd)]; salts = [salts; metData.salts(healyStartInd:healyIncrement:healyEndInd)];
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
SAs = [SAs; nanmean(profiles.SA(1:3, profNums))']; salts = [salts; nanmean(profiles.salts(1:3, profNums))'];
platform = [platform; 2 .* ones(size(profiles.lons(profNums)))];

%Wave Glider data
% wvdata = loadWaveglider;
[SA, ~, ~] = gsw_SA_Sstar_from_SP(wvdata.salts(:, 3), gsw_p_from_z(-wvdata.depths(:,  3), wvdata.lats), wvdata.lons, wvdata.lats);
CT = gsw_CT_from_t(SA, wvdata.temps(:,  3), gsw_p_from_z(-wvdata.depths(:,  3), wvdata.lats));

%Add Wave Glider data to the arrays holding data from all platforms
lons = [lons; wvdata.lons]; lats = [lats; wvdata.lats];
times = [times; wvdata.times];
salts = [salts; wvdata.salts(:, 2)]; %At surface
SAs = [SAs; SA]; %at 9 m so not used
platform = [platform; 3 .* ones(size(wvdata.lons))];

%Load SMOS images
files = dir([rootPath, 'data/SMOS_SSS/*7days.nc']);
if isempty(files_daily)
     disp(['SMOS sea surface salinity data not found. Download weekly daily files for 15 September 2018, 21 September 2018, 30 September 2018, and 7 October 2018 from', ...
            newline, 'https://www.seanoe.org/data/00607/71909/', ...
            newline, 'to the directory ~/meltwaterAdvection/data/SMOS_SSS/'])        
    return
end

curFile = files(1).name;
smosLon = ncread(curFile, 'longitude');
smosLat = ncread(curFile, 'latitude'); %In descending order from 90 to -90

%Iterate through all SMOS images, average SMOS data near each observation
smosSSS_nearObservations = nan .* ones(size(lons));
for timePeriod = 1:length(files)
    curFile = files(timePeriod).name;
    sss = ncread(curFile, 'smos_sss');
    startTime = datenum([curFile(24:27), ',', curFile(29:30), ',' curFile(32:33)]); %Parse string to extract time from file name
    endTime = startTime + 7; %Weekly images
    
    %Iterate through all observations. If in the time range of the current
    %image, average nearby SMOS SSS
    for i = 1:length(lons)
        if times(i) < (startTime-timeThreshold) | times(i) > (endTime + timeThreshold)
            continue; %This observation is not in the time range of the current image. 
        end
        
        %This is the geographic range around each point to search for 
        %SMOS data to average. Box is roughly 15 km
        minlon = lons(i) - 0.25; maxlon = lons(i) + 0.25;  %At 75 N, 0.15 degrees latitude ~= 16 km, 0.5 degrees longitude ~= 14.5 km
        minlat = lats(i) - 0.075; maxlat = lats(i) + 0.075; 
        pt1 = [minlon, minlat]; pt2 = [maxlon, minlat]; pt3 = [maxlon, maxlat]; pt4 = [minlon, maxlat];
        pts = [pt1(1), pt2(1), pt3(1), pt4(1); pt1(2), pt2(2), pt3(2), pt4(2)];

        %Average nearby data
        [inBoxMask, on] = inpolygon(smosLon, smosLat, pts(1, :), pts(2, :));
        smosSSS_nearObservations(i) = nanmean(sss(inBoxMask));
    end
end 

%% Begin plotting
timePeriod = 0; %If timePeriod ~= 1, 2, 3, then will plot data from the entire time on one plot
startTime = datenum('sept 19 2018');
endTime = datenum('oct 15 2018'); 
limits = [23, 30]; %Axis limits for plotting
if timePeriod == 1
    startTime = datenum('sept 19 2018'); 
    endTime = datenum('sept 24 2018');
elseif timePeriod == 2
    startTime = datenum('sept 30 2018'); 
    endTime = datenum('oct 7 2018'); 
elseif timePeriod == 3
    startTime = datenum('oct 8 2018'); 
    endTime = datenum('oct 15 2018'); 
end

inTimeMask = zeros(size(times));
inTimeMask(times >= startTime & times < endTime) = 1;

pedgs = 23:0.2:30;
[X, Y] = meshgrid(pedgs, pedgs);

%% Data from all platforms on one plot

%Number of observations from each platform - use this to decide how to
%subsample
nHealy = length(find(platform == 1));
nProfile = length(find(platform == 2));
nWaveglider = length(find(platform == 3));

%Subset the Healy data
healyInc = round(nHealy / nProfile);
healySubsample = zeros(size(lons));
healySubsample(1:healyInc:nHealy) = 1; healySubsample(nHealy:end) = 1;

%Subset the Wave Glider data
wavegliderInc = 5;
wavegliderSubsample = zeros(size(lons));
wavegliderSubsample(1:(length(lons) - nWaveglider)) = 1; wavegliderSubsample((length(lons) - nWaveglider):wavegliderInc:end) = 1;

%Make a mask identifying which observations to use in the comparison
toPlotMask = zeros(size(lons));
toPlotMask(inTimeMask == 1 & healySubsample == 1 & wavegliderSubsample == 1) = 1;

%Standard deviation of the differences - discussed in text of the paper 
nanstd(salts(toPlotMask == 1) - smosSSS_nearObservations(toPlotMask == 1));

%Make the 2d histogram comparing the in situ observations and SMOS SSS 
binCounts = hist3([salts(toPlotMask == 1), smosSSS_nearObservations(toPlotMask == 1)], 'edges', {pedgs pedgs});
binCounts = binCounts'; %Need to transpose - see "Plot Histogram with Intensity Map" example in hist3 documentation
binCounts(binCounts == 0) = nan;

%Begin plotting
figure(4), set(gcf, 'color', 'w')
subplot(1, 2, 2)
g = pcolor(X, Y, binCounts); set(g, 'linestyle', 'none')
hold on
colormap(gca, brewermap(20, 'greys'))
% title('Data from all platforms', 'fontsize', 12)
xlim(limits); ylim(limits);
x = limits(1):.1:limits(2);
y = x; plot(x, y, 'color', 0.7 .* [1 1 1], 'LineStyle', '-')
grid on; box on
cb = colorbar;
ylabel(cb, 'Number of observations', 'fontsize', 14)
pbaspect([1 1 1])
ylabel('SMOS sea surface salinity (psu)', 'fontsize', 14)
xlabel('Observed sea surface salinity (psu)', 'fontsize', 14)
set(gca, 'fontsize', 12)

if saveFigs
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end

%% Make platform specific comparisons, if there are some observations in the time period
if plotPlatformComparison
    figure(1); set(gcf, 'color', 'w', 'pos', [188 563 1402 385])

    for curPlatform = 1:3

        if sum(platform == curPlatform & inTimeMask == 1) > 0 %There is some data to plot

            switch curPlatform
                case 1
                    cmap = 'Blues';
                    titleText = 'Healy underway data';
                case 2
                    cmap = 'Reds';
                    titleText = 'Seaglider and uCTD data';
                case 3
                    cmap = 'Greens';
                    titleText = 'Wave Glider data';
            end

            subplot(1, 3, curPlatform); hold on
            binCounts = hist3([salts(platform == curPlatform & inTimeMask == 1), smosSSS_nearObservations(platform == curPlatform & inTimeMask == 1)], 'edges', {pedgs pedgs});
            binCounts = binCounts'; %Need to transpose - see "Plot Histogram with Intensity Map" example in hist3 documentation
            binCounts(binCounts == 0) = nan;
            g = pcolor(X, Y, binCounts); set(g, 'linestyle', 'none')
            colormap(gca, brewermap(20, cmap))
            title(titleText, 'fontsize', 12)

            % Note - compare to the hist3 implmentation below to make sure the pcolor plot is correct 
            % hist3([salts(platform == 1), smosSSS_nearObservations(platform == 1)], 'edges', {edgs edgs}, 'CdataMode','auto'); view(2)

        end
    end

    %Format plots
    for splot = 1:3
        subplot(1, 3, splot)
        hold on
        xlim(limits); ylim(limits);
        x = limits(1):.1:limits(2);
        y = x; plot(x, y, 'color', 0.7 .* [1 1 1], 'LineStyle', '-') %Add 1:1 line
        grid on; box on
        cb = colorbar;
        ylabel(cb, 'Number of observations', 'fontsize', 12)
        pbaspect([1 1 1])
        ylabel('SMOS SSS', 'fontsize', 12)
        xlabel('Observed SSS', 'fontsize', 12)
        delete(legend)
    end
end
