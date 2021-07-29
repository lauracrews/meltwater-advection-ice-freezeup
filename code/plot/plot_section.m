% Plots temperature and salinity sections from Seagliders and uCTD. 
% Selects profiles taken by sg198 and sg199 for making composite sections 
clearvars -except rootPath AMSR2 profiles wvdata metData transect
close all

saveFigs = true;
saveDir = [rootPath, 'figures/fig5/'];
saveName = ['section', num2str(transect)];

% transect = 1; %Currently, transect is defined outside of this function, in run_meltwaterAdvection.m

% profiles = loadProfiles;
defineSODAconstants;

mldThreshold = 0.25; %Density change threshold relative to surface for plotting mixed layer depth

%Axis limitst for plotting
climits_temp = [-1.5, 1]; %Color limits for temperature
climits_salt = [25, 29];%34]; %color limits for salinity
ylimits_profiles = [0, 50]; %Depth range for section plots
latlimits = [73.6, 74.45]; %North and south limits for section plots .4

%Set up map projection for the entire SODA-A to SODA-B area - used for
%making maps, but NOT for identifying profiles
minlon_proj = -150; maxlon_proj = -143; minlat_proj = 72; maxlat_proj = 75.5;
m_proj('lambert', 'lon', [minlon_proj maxlon_proj], 'lat', [minlat_proj maxlat_proj]);

%Lat/lon limits used for identifying profiles
minlon = -152; maxlon = -142; minlat = 72; maxlat = 75.7;

%Identify profiles in the geographic region
inRegionMask = zeros(size(profiles.times));    
inRegionMask(profiles.lats >= minlat & profiles.lats <= maxlat ...
    & profiles.lons >= minlon & profiles.lons <= maxlon) = 1;

%Empty arrays to hold data along the transect 
temps = []; salts = []; sigthes = []; mlds = []; lats = []; lons = []; times = [];
for gliderNum = 1:2 %Iterate through Seagliders 198 and 199 
    switch gliderNum
        case 1 %Northern part of transect
            curGlider = 'sg198';
            if transect == 1 %No meltwater
                startTime = datenum('sept 22 2018 06:00:00');
                endTime = datenum('sept 26 2018 12:00:00');
            elseif transect == 2 %Meltwater present - crossing the front            
                startTime = datenum('oct 1 2018 00:06:00');
                endTime = datenum('oct 9 2018 00:00:00');                    
            elseif transect == 3 %Meltwater fills transect
                startTime = datenum('oct 5 2018 21:00:00');
                endTime = datenum('oct 11 2018 01:00:00');
            elseif transect == 4 %After meltwater
                startTime = datenum('oct 11 2018 01:00:00');
                endTime = datenum('oct 15 2018 00:00:00');
            end
            
        case 2 %Southern part of transect
            curGlider = 'sg199';
            if transect == 1 %No meltwater
                startTime = datenum('sept 22 2018 00:00:00');
                endTime = datenum('sept 26 2018 20:00:00');
            elseif transect == 2 %Meltwater present - crossing the front        
                startTime = datenum('sept 30 2018 06:00:00');
                endTime = datenum('oct 5 2018 20:00:00');
            elseif transect == 3  %Meltwater fills transect 
                startTime = datenum('oct 5 2018 21:00:00');
                endTime = datenum('oct 11 2018 01:00:00');
            elseif transect == 4 %After meltwater
                startTime = datenum('oct 11 2018 01:00:00');
                endTime = datenum('oct 15 2018 00:00:00');
            end            
    end
    
    %Identify profiles in the time period collected by the correct glider
    inTimeMask = zeros(size(profiles.times));
    inTimeMask(profiles.times >= startTime & profiles.times < endTime) = 1;

    vehicleMask = strcmp(profiles.dataset, curGlider);   
    
    profNums = find(profiles.qualFlag == 1 & inTimeMask == 1 & inRegionMask == 1 & vehicleMask == 1);
    
    %Calculate mixed layer depth of each profile
    for i = 1:length(profiles.lats)
        mldind = find(profiles.sigthes(:, i) - profiles.sigthes(1, i) >= mldThreshold, 1, 'first');
        profiles.mld(i) = profiles.z(mldind);
    end

    %Add data to the growing arrays for the whole transect
    temps = [temps, profiles.CT(:, profNums)]; salts = [salts, profiles.SA(:, profNums)]; sigthes = [sigthes, profiles.sigthes(:, profNums)]; 
    mlds = [mlds, profiles.mld(profNums)];
    lats = [lats, profiles.lats(profNums)']; lons = [lons, profiles.lons(profNums)']; times = [times, profiles.times(profNums)'];
    
    %Key to determine which glider 
    if gliderNum == 1
        gliderIDs = ones(size(profNums'));
    else
        gliderIDs = [gliderIDs, 2 .* ones(size(profNums'))];
    end
end

% %Update arrays to be sorted by latitude
[~, idx] = sort(lats);
temps = temps(:, idx); salts = salts(:, idx); sigthes = sigthes(:, idx); 
lats = lats(idx); lons = lons(idx); times = times(idx);
mlds = mlds(idx); gliderIDs = gliderIDs(idx);

%% Identify the northernmost and southernmost profiles taken by each glider, 
%to be used for labelling the dates/gliders
inds(1) = find(lats == min(lats(gliderIDs == 1 & lats >= latlimits(1) & lats <= latlimits(2))));
inds(2) = find(lats == max(lats(gliderIDs == 1 & lats >= latlimits(1) & lats <= latlimits(2)))); 
inds(3) = find(lats == min(lats(gliderIDs == 2 & lats >= latlimits(1) & lats <= latlimits(2))));
inds(4) = find(lats == max(lats(gliderIDs == 2 & lats >= latlimits(1) & lats <= latlimits(2))));

%% Make maps showing the locations of all profiles included in the section

figure(1); set(gcf, 'color', 'w', 'pos', [144 227 1449 728])

%Add MODIS SST imagery to the background
if transect == 1
    curFile = 'T20182572018264.L3m_8D_NSST_sst_4km.nc'; 
    iceDay = datenum('Sept 25 2018');
elseif transect == 2
    curFile = 'T20182732018280.L3m_8D_NSST_sst_4km.nc';
    iceDay = datenum('Oct 4 2018');
end
[modisLon, modisLat, sst, ~, ~, modisTitle] = loadMODISsst(curFile, minlon_proj, maxlon_proj, minlat_proj, maxlat_proj);

%First axis is for sea surface temperature
ax1 = axes; subplot(1, 4, 1, ax1)
h1 = m_pcolor(modisLon, modisLat, sst); %MODIS background sst
hold on 

%Add in sity obsevations to the map
m_scatter(lons(lats >= latlimits(1) & lats <= latlimits(2)), lats(lats >= latlimits(1) & lats <= latlimits(2)), 80, 'w', 'filled');

m_scatter(lons(lats >= latlimits(1) & lats <= latlimits(2)), lats(lats >= latlimits(1) & lats <= latlimits(2)), 60,...
    nanmean(temps(1:5, (lats >= latlimits(1) & lats <= latlimits(2)))), 'filled');

cmocean('thermal'); caxis(climits_temp)  
m_grid('xtick', [-149:2:-143], 'linestyle', 'none', 'fontsize', 12) %Need to m_grid on both sets of axes 
ps = get(gca, 'pos');
cb = colorbar;
set(cb, 'fontsize', 12, 'location', 'southoutside')
ylabel(cb, ['Conservative temperature ', ' (', sprintf(char(176)), 'C)'], 'fontsize', 16)%
m_gshhs_i('patch',0.6*[1 1 1]); %Draw coastline

obsTitle = ['Observations ' datestr(min(times(inds)), 'mmm dd'), ' to ', datestr(max(times(inds)), 'mmm dd')];
iceTitle = ['AMSR2 sea ice ' datestr(iceDay, 'mmm dd')];
title([modisTitle, sprintf('\n'), obsTitle, sprintf('\n'), iceTitle], 'fontsize',14, 'fontweight', 'bold') 

% Add another axes for the ice concentration. New axis is necessary to allow different colormap
ax2 = axes;
subplot(1, 4, 1, ax2)  
[~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
curIce = nanmean(AMSR2.SIC(:, :, curIceInd), 3);
curIce(curIce <= 0) = nan;
h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
set(h2, 'facealpha', .8)
cmocean('ice')
hold on
m_grid('xticklabels', [], 'yticklabels', [], 'linestyle', 'none')
cb2 = colorbar;
cbps = get(cb, 'pos'); cbps(2) = cbps(2) - 0.13;
set(cb2, 'fontsize', 12, 'location', 'southoutside')
set(cb2, 'position', cbps)
ylabel(cb2, [' AMSR2 sea ice concentration'], 'fontsize', 12)%  
set(ax2, 'pos', ps)
set(ax1, 'pos', ps)

m_scatter(moorings([1, 3:4], 2), moorings([1, 3:4], 1), 250, 'k', 'p', 'filled')
m_scatter(moorings([1, 3:4], 2), moorings([1, 3:4], 1), 200, 'w', 'p', 'filled')

%Save map figure
if saveFigs  
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName, '_map'],'-dpng')
    saveas(gcf, [saveDir, saveName, '_map.fig'])
end

%% Make section plots

%Make grids of latitude and depth, for plotting sections
[LATS, DEPTHS] = meshgrid(lats, profiles.z);

figure(2); set(gcf, 'pos', [244 366 1135 582], 'color', 'w')

subplot(2, 1, 1) %Plot temperature section
pcolor(LATS, DEPTHS, temps); shading interp
hold on
cmocean('thermal'); caxis(climits_temp)   
ps = get(gca, 'pos');
cb = colorbar;
ylabel(cb, ['Conservative temperature (', sprintf(char(176)), 'C)'], 'fontsize', 14) 
set(gca, 'fontsize', 12, 'pos', ps)
title(['Profiles taken between ', datestr(min(times(inds)), 'mmm dd'), ' and ', datestr(max(times(inds)), 'mmm dd yyyy')], 'fontsize', 14)

subplot(2, 1, 2) %Plot salinity section
pcolor(LATS, DEPTHS, salts); shading interp
hold on
cmocean('haline'); caxis(climits_salt)
ps = get(gca, 'pos');
cb = colorbar;
ylabel(cb, 'Absolute salinity (g/kg)', 'fontsize', 14) 
xlabel('Latitude', 'fontsize', 14) 
set(gca, 'fontsize', 12, 'pos', ps)

%Label the part of the transect completed by each glider 
text(lats(inds), -3.* ones(size(inds)), datestr(times(inds), 'mmm dd'), 'fontsize', 12, 'rotation', 90)
text(mean([lats(inds(1)), lats(inds(2))]) -.04, -8, 'Seaglider 198', 'fontsize', 14) %Northern glider  
text(mean([lats(inds(3)), lats(inds(4))]) -.04, -8, 'Seaglider 199', 'fontsize', 14) %Northern glider  

%Format section plots
for splot = 1:2
    subplot(2, 1, splot)
    [C, h] = contour(LATS, DEPTHS, sigthes, [1020:.25:1028], 'color', 0.6 .* [1 1 1]); %Add density contours
    clabel(C, h, 'color', 0.6 .* [1 1 1]) %Label density contours
    plot(lats, mlds, 'w', 'linewidth', 1, 'linestyle', '--') %Add mixed layer depth
    set(gca, 'YDir', 'reverse')
    ylim(ylimits_profiles); xlim(latlimits)
    ylabel('Depth (m)', 'fontsize', 14)

    for j = 1:length(lons) %Add lines showing the location of each profile
        line([lats(j), lats(j)], ylimits_profiles, 'color', 'k', 'linestyle', ':')
    end
end

%Save figures
if saveFigs    
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end