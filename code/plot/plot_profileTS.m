%Makes temperature-salinity plots for profiles within a given geographic
%and time range. Also makes a map showing where the profiles were taken.
close all; 
clearvars -except rootPath AMSR2 profiles wvdata metData 
saveFigs = false;

fntsz = 16;
saveDir = [rootPath, 'figures/fig6/'];
saveName = 'meltwaterTS';

%Load Seaglider and uCTD profiles
% profiles = loadProfiles;

%Axis limits for plotting
clim_temp = [-1.5, 1];
xlimits_temp = [-1.9, 2];
xlimits_salt = [25, 35];
[moorings, colors] = defineSODAconstants;

%% Map SST to overlay profile locations. 
%Need to make map before iterating through and adding profiles

%Set up map projection
minlon = -150; maxlon = -143; minlat = 72; maxlat = 75.5;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

%Will add sea ice to map
iceDay = datenum('Oct 4 2018');
iceTitle = ['AMSR2 sea ice ' datestr(iceDay, 'dd mmm')];

%MODIS sea surface temperature image to be the background of the map
curFile = 'T20182732018280.L3m_8D_NSST_sst_4km.nc';
[modisLon, modisLat, sst, ~, ~, modisTitle] = loadMODISsst(curFile, minlon, maxlon, minlat, maxlat);

%First axis is for the MODIS SST image
% %figure(2); set(gcf, 'color', 'w', 'pos', [144 227 1449 728])

ax1 = axes; subplot(1, 4, 1, ax1)
h1 = m_pcolor(modisLon, modisLat, sst);
hold on  
set(h1, 'facealpha', .9)
caxis(clim_temp); cmocean('thermal')
freezeColors(h1)
ps = get(gca, 'pos');
cb = colorbar;
m_grid('xtick', [-149:2:-143], 'fontsize', fntsz)
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
m_grid('xticklabels', [], 'yticklabels', [], 'linestyle', 'none') 
cb2 = colorbar;
cbps = get(cb, 'pos'); cbps(2) = cbps(2) - 0.22;%0.13;
set(cb2, 'location', 'southoutside') %Both of these cb positioning lines are necessary!
set(cb2, 'fontsize', 14, 'position', [0.13    0.1   0.1562    0.0340])%cbps)

set(cb, 'location', 'southoutside'); 
set(cb, 'fontsize', 14, 'pos', [0.13    0.21    0.1562    0.034])

ylabel(cb, ['Conservative temperature ', ' (', sprintf(char(176)), 'C)'])%
ylabel(cb2, ['AMSR2 sea ice concentration'])%  
set(ax2, 'pos', ps)
set(ax1, 'pos', ps)
m_scatter(moorings.all(:, 2), moorings.all(:, 1), 250, 'k', 'p', 'filled')
m_scatter(moorings.all(:, 2), moorings.all(:, 1), 200, 'w', 'p', 'filled')

%% Begin making Figure 1 which will be T-S plot
%figure(1); set(gcf, 'Position', [407   284   810   664], 'color', 'w') 
set(gcf, 'Position', [224 281 1188 666], 'color', 'w') 
subplot(1, 5, 3:5)
hold on; box on
tempGrid = -2:0.01:2; saltGrid = 24:0.01:36;
freezePt_TSplot = gsw_CT_freezing(saltGrid, 0);
plot(saltGrid, freezePt_TSplot, 'k')
text(25.2, -1.55, 'Freezing temperature', 'fontsize', 16, 'FontWeight', 'bold')

% Calculate and plot Gade lines (Gade, 1979, Journal of Physical Oceanography)
T1 = freezePt_TSplot;
T0 = -2 .* ones(size(T1));
cw = gsw_cp_t_exact(saltGrid, T1, zeros(size(T1)));
ci = gsw_cp_ice(T1, 0);

T3 = 2 .* ones(size(T1));
S3 = saltGrid;
L = gsw_latentheat_melting(saltGrid, zeros(size(saltGrid)));

slope = (1 ./ S3) .* (T3 - T1 + (T1-T0)*(ci/cw) + L./cw);

T2 = T1; %on the freezing line
S2 = saltGrid;
S3 = S2 + (T3 - T2)./ slope;
for i = 150:30:250
    line([S3(i), S2(i)], [T3(i), T2(i)], 'linewidth', 1, 'color', [1.00,0.41,0.16]);
end


%%
[saltGrid, tempGrid] = meshgrid(saltGrid, tempGrid);

sigtheGrid = gsw_rho(saltGrid, tempGrid, 0);

%Mask out density values below the freezing temperature
for colNum = 1:size(sigtheGrid, 2)
    sigtheGrid(tempGrid(:, colNum) < freezePt_TSplot(colNum), colNum) = nan;
end

[C, h] = contour(saltGrid, tempGrid, sigtheGrid, 'k', 'LineWidth', 1, 'LineStyle', ':');
clabel(C, h, [1021, 1023, 1025, 1027, 1029], 'fontsize', 12)



%Bold the 23.2 and 25.2 isopycnals above 0C and add the 0C line, per the PSW definition used in MacKinnon et al. (2021)
pswTempGrid = tempGrid; pswTempGrid(tempGrid<0) = nan;
pswSaltGrid = saltGrid; pswSaltGrid(tempGrid<0) = nan;
pswSigtheGrid = sigtheGrid; pswSigtheGrid(tempGrid<0) = nan;
contour(pswSaltGrid, pswTempGrid, pswSigtheGrid, [1023.3, 1025.2], 'k', 'LineWidth', 2, 'LineStyle', '-');
line([29.17 31.56], [0 0], 'color', 'k', 'LineWidth', 2, 'LineStyle', '-');

%Label water masses
text(25.2, -1, 'MW', 'fontsize', 14, 'FontWeight', 'bold')
text(29.5, 1.7, 'PSW', 'fontsize', 14, 'FontWeight', 'bold')
text(32.8, -1.6, 'PWW', 'fontsize', 14, 'FontWeight', 'bold')
text(34, -0.2, 'AW', 'fontsize', 14, 'FontWeight', 'bold')

xlim([xlimits_salt])
xlabel(['Absolute salinity (g/kg)'], 'fontsize', 16)%
ylabel(['Conservative temperature ', ' (', sprintf(char(176)), 'C)'], 'fontsize', 16)%
set(gca, 'fontsize', 14)
%%
plotCount = 1;
for group = [3, 1, 2] %Iterate through combinations of lat/lon and time bounds
    switch group
        case 1 %Before meltwater
            minlon = -148; maxlon = -146; minlat = 73 + 25/60; maxlat = 74.35; 
            startTime = datenum('sept 22 2018'); endTime = datenum('oct 3 2018');
            col = colors.blue;
        case 2 %Meltwater present
            minlon = -148; maxlon = -146; minlat = 73 + 25/60; maxlat = 74.35;
            startTime = datenum('oct 5 2018'); endTime = datenum('oct 11 2018');
            col = colors.purple;
        case 3 %Remnant ice area
            minlon = -149.5; maxlon = -147; minlat = 72.75; maxlat = 73.25;
            startTime = datenum('sept 26, 2018'); endTime = datenum('sept 26, 2018') + .95;
            col = 'k';
    end
    
    %Identify the profiles in this geographic region in this time period    
    profNums = identifyProfileNumbers(profiles, minlon, minlat, maxlon, maxlat, startTime, endTime);
     
    %Plot profiles on T-S plot
    %figure(1);
    subplot(1, 5, 3:5)
    for i = 1:length(profNums)
        h3(plotCount) = scatter(profiles.SA(:, profNums(i)), profiles.CT(:, profNums(i)), 4, col);
    end
    
    % Plot profiles on map
    %figure(2); 
    
    subplot(1, 4, 1, ax2)
    h(plotCount) = m_scatter(profiles.lons(profNums), profiles.lats(profNums), 20, col);
    hold on 
    if group ~= 3
        lgdText{plotCount} = ['Profiles taken ', datestr(datestr(startTime), 'dd mmm'), ' to ', datestr(datestr(endTime), 'dd mmm')];
    else
        lgdText{plotCount} = ['Profiles taken ', datestr(datestr(startTime), 'dd mmm')];
    end
   
    plotCount = plotCount + 1;
end


%% Format and save plots
%figure(1)
box on

lgd = legend(h(1:3), lgdText{1}, lgdText{2}, lgdText{3});
set(lgd, 'fontsize', 14, 'pos', [0.1069    0.8211    0.2012    0.0818])

if saveFigs 
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end

%figure(2)
if saveFigs    
    print([saveDir, saveName, '_map'],'-dpng')
    saveas(gcf, [saveDir, saveName, '_map.fig'])
end
