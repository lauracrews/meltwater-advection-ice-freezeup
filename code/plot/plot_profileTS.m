%Makes temperature-salinity plots for profiles within a given geographic
%and time range. Also makes a map showing where the profiles were taken.
close all; 
clearvars -except rootPath AMSR2 profiles wvdata metData 

saveFigs = true;
saveDir = [rootPath, 'figures/fig6/'];
saveName = 'meltwaterTS';

%Load Seaglider and uCTD profiles
% profiles = loadProfiles;

%Axis limits for plotting
clim_temp = [-1.5, 1];
xlimits_temp = [-1.9, 2];
xlimits_salt = [25, 34];
defineSODAconstants;

%% Map SST to overlay profile locations. 
%Need to make map before iterating through and adding profiles

%Set up map projection
minlon = -150; maxlon = -143; minlat = 72; maxlat = 75.5;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

%Will add sea ice to map
load AMSR2_2018.mat
iceDay = datenum('Oct 4 2018');
iceTitle = ['AMSR2 sea ice ' datestr(iceDay, 'mmm dd')];

%MODIS sea surface temperature image to be the background of the map
curFile = 'T20182732018280.L3m_8D_NSST_sst_4km.nc';
[modisLon, modisLat, sst, ~, ~, modisTitle] = loadMODISsst(curFile, minlon, maxlon, minlat, maxlat);

%First axis is for the MODIS SST image
figure(2); set(gcf, 'color', 'w', 'pos', [144 227 1449 728])
ax1 = axes; subplot(1, 4, 1, ax1)
h1 = m_pcolor(modisLon, modisLat, sst);
hold on  
set(h1, 'facealpha', .9)
caxis(clim_temp); cmocean('thermal')
freezeColors(h1)
ps = get(gca, 'pos');
cb = colorbar;
m_grid('fontsize', 12)
set(cb, 'fontsize', 12, 'location', 'southoutside')
ylabel(cb, ['Conservative temperature ', ' (', sprintf(char(176)), 'C)'], 'fontsize', 16)%
m_gshhs_i('patch',0.6*[1 1 1]); %Draw coastline
title([modisTitle, sprintf('\n'), iceTitle], 'fontsize',14, 'fontweight', 'bold') 

% Add another axes for the ice concentration. New axis is necessary to allow different colormap
ax2 = axes; subplot(1, 4, 1, ax2)
[~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
curIce = nanmean(AMSR2.SIC(:, :, curIceInd), 3);
curIce(curIce <= 0) = nan;
h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
set(h2, 'facealpha', .8)
cmocean('ice')
freezeColors(h2)
hold on
m_grid('linestyle', 'none', 'fontsize', 12) %Need to m_grid on both sets of axes 
cb2 = colorbar;
cbps = get(cb, 'pos'); cbps(2) = cbps(2) - 0.13;
set(cb2, 'fontsize', 12, 'location', 'southoutside') %Both of these cb positioning lines are necessary!
set(cb2, 'position', cbps, 'fontsize', 12)
ylabel(cb2, ['AMSR2 sea ice concentration'], 'fontsize', 16)%  
set(ax2, 'pos', ps)
set(ax1, 'pos', ps)
m_scatter(moorings(:, 2), moorings(:, 1), 250, 'k', 'p', 'filled')
m_scatter(moorings(:, 2), moorings(:, 1), 200, 'w', 'p', 'filled')

%% Begin making Figure 1 which will be T-S plot
figure(1); set(gcf, 'Position', [407   284   810   664], 'color', 'w') 
hold on
tempGrid = -2:0.01:2; saltGrid = 24:0.01:34;
freezePt_TSplot = gsw_CT_freezing(saltGrid, 0);
plot(saltGrid, freezePt_TSplot, 'k')
text(25.5, -1.75, 'Freezing temperature', 'fontsize', 16, 'FontWeight', 'bold')
[saltGrid, tempGrid] = meshgrid(saltGrid, tempGrid);
sigtheGrid = gsw_rho(saltGrid, tempGrid, 0);
[C, h] = contour(saltGrid, tempGrid, sigtheGrid, 'k', 'LineWidth', 1, 'LineStyle', ':');
clabel(C, h, [1021, 1023, 1025, 1027, 1029])
xlabel(['Absolute salinity (g/kg)'], 'fontsize', 14)%
ylabel(['Conservative temperature ', ' (', sprintf(char(176)), 'C)'], 'fontsize', 14)%
set(gca, 'fontsize', 12)

plotCount = 1;
for group = [1:2, 8] %Iterate through combinations of lat/lon and time bounds
    switch group
        case 1 %South, before meltwater
            minlon = -148; maxlon = -146; minlat = 73 + 25/60; maxlat = 74.35; 
            startTime = datenum('sept 22 2018'); endTime = datenum('oct 3 2018');
            col = blue;
        case 2 %South, meltwater present
            minlon = -148; maxlon = -146; minlat = 73 + 25/60; maxlat = 74.35;
            startTime = datenum('oct 5 2018'); endTime = datenum('oct 11 2018');
            col = purple;
        case 3 %South, after meltwater
            minlon = -148; maxlon = -146; minlat = 73 + 25/60; maxlat = 74.35;
            startTime = datenum('oct 11 2018'); endTime = datenum('oct 15 2018'); 
            col = ltblue;
        case 4 %Mid latitudes (indludes eddy)
            minlon = -147; maxlon = -145.5; minlat = 74.35; maxlat = 74.67;
            startTime = datenum('sept 21 2018'); endTime = datenum('oct 15 2018'); %endTime = datenum('sept 27 2018');
            col = orange;
        case 5 %Mid latitudes
            minlon = -147; maxlon = -145.5; minlat = 74.35; maxlat = 74.67;
            startTime = datenum('sept 27 2018'); endTime = datenum('oct 6 2018');
            col = orange;
        case 6 %Northern end
            minlon = -147; maxlon = -145; minlat = 74.67; maxlat = 75.1;
            startTime = datenum('sept 22 2018'); endTime = datenum('oct 15 2018'); %endTime = datenum('sept 30 2018'); 
            col = red;
         case 7 %Northern end
            minlon = -147; maxlon = -145; minlat = 74.67; maxlat = 75.1;
            startTime = datenum('sept 30 2018'); endTime = datenum('oct 8 2018');
            col = orange;
        case 8 %Remnant ice area
            minlon = -149.5; maxlon = -147; minlat = 72.75; maxlat = 73.25;
            startTime = datenum('sept 26, 2018'); endTime = datenum('sept 26, 2018') + .95;
            col = 'k';
        case 9
            minlon = -147; maxlon = -145; minlat = 74.67; maxlat = 75.1;
            startTime = datenum('oct 8 2018'); endTime = datenum('oct 15 2018'); %endTime = datenum('sept 30 2018'); 
            col = red;
        case 10
            minlon = -147; maxlon = -145; minlat = 75.1; maxlat = 75.5;
            startTime = datenum('oct 11 2018'); endTime = datenum('oct 12 2018'); %endTime = datenum('sept 30 2018'); 
            col = blue;
    end
    
    %Mask of profiles in the current geographic bounds
    inRegionMask = zeros(size(profiles.times));    
    inRegionMask(profiles.lats >= minlat & profiles.lats <= maxlat ...
        & profiles.lons >= minlon & profiles.lons <= maxlon) = 1;

    %Mask of profiles in the current time bounds
    inTimeMask = zeros(size(profiles.times));
    inTimeMask(profiles.times >= startTime & profiles.times < endTime) = 1;
     
    %Identify the profiles in this geographic region in this time period
    if group ~=3
        profNums = find(profiles.qualFlag == 1 & inTimeMask == 1 & inRegionMask == 1); 
    else
        vehicleMask = strcmp(profiles.dataset, 'uCTD');
        profNums = find(profiles.qualFlag == 1 & inTimeMask == 1 & inRegionMask == 1 & vehicleMask == 1); 
    end
    
    %Plot profiles on T-S plot
    figure(1);
    for i = 1:length(profNums)
        h3(plotCount) = scatter(profiles.SA(:, profNums(i)), profiles.CT(:, profNums(i)), 4, col);
    end
    
    % Plot profiles on map
    figure(2); 
    h(plotCount) = m_scatter(profiles.lons(profNums), profiles.lats(profNums), 20, col);
    hold on 
    lgdText{plotCount} = ['Profiles taken ', datestr(datestr(startTime), 'mmm dd'), ' to ', datestr(datestr(endTime), 'mmm dd')];
   
    plotCount = plotCount + 1;
end
%% Format and save plots
figure(1)
xlim([25, 33])
box on
if saveFigs 
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end

figure(2)
lgd = legend(h(1:3), lgdText{1}, lgdText{2}, lgdText{3});
set(lgd, 'fontsize', 12)
if saveFigs    
    print([saveDir, saveName, '_map'],'-dpng')
    saveas(gcf, [saveDir, saveName, '_map.fig'])
end
