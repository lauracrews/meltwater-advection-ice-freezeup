% Plots temperature and salinity sections from Seagliders and uCTD. 
% Selects profiles taken by sg198 and sg199 for making composite sections 
clearvars -except rootPath AMSR2 profiles wvdata metData transect saveFigs
% close all

% transect = 1; %Currently, transect is defined outside of this function, in run_meltwaterAdvection.m

saveFigs = false;
saveDir = [rootPath, 'figures/fig5/'];
saveName = ['section', num2str(transect)];
fntsz = 16; %font size for entire figure. Nice to be able to update here - want big font for paper but smaller is better for presentations

makeCompositeFigure = true; %If all panels will be combined into one figure, skip some labels and colorbars to reduce text redundancy
labelSampler = true; %Label which glider collected each part of the transect

[moorings, ~] = defineSODAconstants;

mldThreshold = 0.25; %Density change threshold relative to surface for plotting mixed layer depth

%Axis limits for plotting
climits_temp = [-1.5, 1]; %Color limits for temperature
climits_salt = [25, 29];%34]; %color limits for salinity
ylimits_profiles = [0, 50]; %Depth range for section plots
latlimits = [73.615, 74.45]; %North and south limits for section plots .4

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
latlimits = [min(lats), max(lats)]; %North and south limits for section plots .4


if labelSampler
    inds(1) = find(lats == min(lats(gliderIDs == 1 & lats >= latlimits(1) & lats <= latlimits(2))));
    inds(2) = find(lats == max(lats(gliderIDs == 1 & lats >= latlimits(1) & lats <= latlimits(2)))); 
    inds(3) = find(lats == min(lats(gliderIDs == 2 & lats >= latlimits(1) & lats <= latlimits(2))));
    inds(4) = find(lats == max(lats(gliderIDs == 2 & lats >= latlimits(1) & lats <= latlimits(2))));
end

%% Make maps showing the locations of all profiles included in the section

alphabet = ('a':'z').';
chars = num2cell(alphabet);
chars = chars.';

figure(1); 
set(gcf, 'color', 'w', 'pos', [9          91        1672         856])

%Add MODIS SST imagery to the background
if transect == 1
    curFile = 'T20182572018264.L3m_8D_NSST_sst_4km.nc'; 
    iceDay = datenum('Sept 18 2018');
elseif transect == 2
    curFile = 'T20182732018280.L3m_8D_NSST_sst_4km.nc';
    iceDay = datenum('Oct 4 2018');
elseif transect == 3;
    iceDay = datenum('Oct 8 2018');
    modisTitle = '';
end

%First axis is for sea surface temperature
ax1 = axes; subplot(2, 5, 5*transect-3, ax1)
hold on 

%Add MODIS SST map, if available
if transect < 3
    [modisLon, modisLat, sst, ~, ~, modisTitle] = loadMODISsst(curFile, minlon_proj, maxlon_proj, minlat_proj, maxlat_proj);
    h1 = m_pcolor(modisLon, modisLat, sst); %MODIS background sst
end

%Add in sity obsevations to the map
m_scatter(lons(lats >= latlimits(1) & lats <= latlimits(2)), lats(lats >= latlimits(1) & lats <= latlimits(2)), 80, 'w', 'filled');

m_scatter(lons(lats >= latlimits(1) & lats <= latlimits(2)), lats(lats >= latlimits(1) & lats <= latlimits(2)), 60,...
    nanmean(temps(1:5, (lats >= latlimits(1) & lats <= latlimits(2)))), 'filled');

cmocean('thermal'); caxis(climits_temp)  
m_grid('xtick', [-149:2:-143], 'linestyle', 'none', 'fontsize', fntsz) %Need to m_grid on both sets of axes 
ps = get(gca, 'pos');

if transect == 1
    cbtemp = colorbar;
    set(cbtemp, 'fontsize', fntsz, 'location', 'southoutside')
    ylabel(cbtemp, ['Conservative temperature ', ' (', sprintf(char(176)), 'C)'], 'fontsize', fntsz)%
end
m_gshhs_i('patch',0.6*[1 1 1]); %Draw coastline

% ['Observations ' datestr(min(times(inds)), 'dd mmm'), ' to ', datestr(max(times(inds)), 'dd mmm')];
obsTitle = ['Observations ', datestr(min(times(lats>=latlimits(1) & lats<=latlimits(2))), 'dd mmm'), ' to ', datestr(max(times(lats>=latlimits(1) & lats<=latlimits(2))), 'dd mmm')];
iceTitle = ['AMSR2 sea ice ' datestr(iceDay, 'dd mmm')];
title([modisTitle, sprintf('\n'), iceTitle], 'fontsize',fntsz-2)%, 'fontweight', 'bold') %obsTitle, sprintf('\n'),

% Add another axes for the ice concentration. New axis is necessary to allow different colormap
ax2 = axes;
subplot(2, 5, 5*transect-3, ax2)  
[~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
curIce = nanmean(AMSR2.SIC(:, :, curIceInd), 3);
curIce(curIce <= 0) = nan;
h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
% set(h2, 'facealpha', .8)
cmocean('ice')
hold on
m_grid('xticklabels', [], 'yticklabels', [], 'linestyle', 'none')

if transect == 1
    cbice = colorbar;
    cbps = get(cbtemp, 'pos'); cbps(2) = cbps(2) - 0.13;
    set(cbice, 'fontsize', fntsz, 'location', 'southoutside')
    set(cbice, 'position', cbps)
    ylabel(cbice, [' AMSR2 sea ice concentration'], 'fontsize', fntsz)%
end
set(ax2, 'pos', ps)
set(ax1, 'pos', ps)

m_scatter(moorings.all([1, 3:4], 2), moorings.all([1, 3:4], 1), 250, 'k', 'p', 'filled')
m_scatter(moorings.all([1, 3:4], 2), moorings.all([1, 3:4], 1), 200, 'w', 'p', 'filled')

if makeCompositeFigure & transect == 1
    set(cbtemp, 'pos', [0.096    0.74    0.1500    0.0320])%[0.34 0.76 0.15 0.032]) 
    set(cbice, 'pos', [0.0960    0.86    0.1500    0.0320])
end

%Subplot letter label
text(0.125,0.95,chars{3*transect - 2},'Units','normalized','color', 'w', 'FontSize',fntsz, 'fontweight', 'bold')

%Save map figure
if saveFigs  
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName, '_map'],'-dpng')
    saveas(gcf, [saveDir, saveName, '_map.fig'])
end

%% Make section plots
n = 30; %Number of distance ticks
ind = find(lats > latlimits(2), 1, 'first') - 1; if isempty(ind); [~, ind] = max(lats); end 
lonpts = linspace(lons(1), lons(ind), n);
latpts = linspace(lats(1), lats(ind), n);
delta = m_lldist(lonpts, latpts); %Displacement along section
dists = cumsum([0; delta])'; %Cumulative distance along section

dists2 = [0:5:max(dists)];
lats2 = interp1(dists, latpts, dists2);

%Make grids of latitude and depth, for plotting sections
[LATS, DEPTHS] = meshgrid(lats, profiles.z);

% figure(2); set(gcf, 'pos', [244 366 1135 582], 'color', 'w')

ax1 = axes;

if transect == 1 
    subplot(4, 5, 5*transect-2:5*transect)
elseif transect == 2
    subplot(4, 5, 5*transect+3:5*transect+5) %Plot temperature section
end

pcolor(LATS, DEPTHS, temps); shading interp
hold on
cmocean('thermal'); caxis(climits_temp)   

if makeCompositeFigure & transect == 1
    for i = 1:length(dists2)
        text(lats2(i), -3, num2str(round(dists2(i))), 'fontsize', fntsz)
    end
    title(['Distance along section (km)', newline], 'fontsize', fntsz, 'FontWeight', 'normal')
end

if ~makeCompositeFigure %Don't need to title the section panels because the maps are already titled
%     title(['Profiles taken between ', datestr(min(times(lats>=latlimits(1) & lats<=latlimits(2))), 'dd mmm'), ' and ', datestr(max(times(lats>=latlimits(1) & lats<=latlimits(2)))), 'dd mmm')], 'fontsize', fntsz)
    title(['Profiles taken between ', datestr(min(times(lats>=latlimits(1) & lats<=latlimits(2))), 'dd mmm'), ' and ', datestr(max(times(lats>=latlimits(1) & lats<=latlimits(2))), 'dd mmm')], 'fontsize', fntsz)

    %Since color limits are constant throughout the entire figure, don't need another colorbar
    ps = get(gca, 'pos');
    cb = colorbar('location', 'eastoutside');
    ylabel(cb, ['Conservative temperature (', sprintf(char(176)), 'C)'], 'fontsize', fntsz) 
    set(gca, 'pos', ps)
else
    set(gca, 'xtick', [])
    text(0.025,0.9,chars{3*transect - 1},'Units','normalized','color', 'w', 'FontSize',fntsz, 'fontweight', 'bold')
end

subplot(4, 5, 10*transect-2:10*transect) %Plot salinity section
pcolor(LATS, DEPTHS, salts); shading interp
hold on
cmocean('haline'); caxis(climits_salt)

if makeCompositeFigure
    text(0.025,0.9,chars{3*transect},'Units','normalized','color', 'w', 'FontSize', fntsz, 'fontweight', 'bold')
end

if makeCompositeFigure & transect == 1   
    set(gca, 'xtick', [])
    ps = get(gca, 'pos');
    if makeCompositeFigure 
        cbsalt = colorbar('location', 'southoutside');
        set(cbsalt, 'pos', [0.096    0.62    0.1500    0.0320])
    else
        cbsalt = colorbar('location', 'eastoutside');
    end
    ylabel(cbsalt, 'Absolute salinity (g/kg)', 'fontsize', fntsz) 
    set(gca, 'pos', ps)
end

if transect == 2
    xlabel(['Latitude (', sprintf(char(176)), 'N)'], 'fontsize', fntsz) 
end
    
%Label the part of the transect completed by each glider 
if labelSampler
    text(lats(inds), -3.* ones(size(inds)), datestr(times(inds), 'dd mmm'), 'fontsize', fntsz-4, 'rotation', 90)
    text(mean([lats(inds(1)), lats(inds(2))]) -.04, -8, 'Seaglider 198', 'fontsize', fntsz-4) %Northern glider  
%     text(mean([lats(inds(3)), lats(inds(4))]) -.04, -8, 'Seaglider 199', 'fontsize', fntsz-4) %Southern glider  
    text(73.92, -8, 'Seaglider 199', 'fontsize', fntsz-4) %Southern glider  
end

%Format section plots
for splot = 1:2
    switch splot
        case 1
            if transect == 1 
                subplot(4, 5, 5*transect-2:5*transect)
            elseif transect == 2
                subplot(4, 5, 5*transect+3:5*transect+5) %Plot temperature section
            end
%             line([73.92, 74], 40 .* [1 1], 'linewidth', 1, 'color', 'k')

        case 2
            subplot(4, 5, 10*transect-2:10*transect) 

    end
    [C, h] = contour(LATS, DEPTHS, sigthes, [1020:.5:1028], 'color', 0.6 .* [1 1 1]); %Add density contours
    clabel(C, h, 'color', 0.6 .* [1 1 1], 'fontsize', fntsz-2) %Label density contours
    plot(lats, mlds, 'w', 'linewidth', 1, 'linestyle', '--') %Add mixed layer depth
    set(gca, 'YDir', 'reverse')
    ylim(ylimits_profiles); xlim(latlimits)
    ylabel('Depth (m)', 'fontsize', fntsz)
    set(gca, 'fontsize', fntsz)

    for j = 1:length(lons) %Add lines showing the location of each profile
        line([lats(j), lats(j)], ylimits_profiles, 'color', 'k', 'linestyle', ':')
    end
end

%Save figures
if saveFigs   
        
%     Lowers the upper panel to squeeze subplots. Seems to work if you
%     run it in the command line once the figure is done, but not if you
%     just include the code (then the figure goes blank)
%     set(gca, 'pos', [0.1300    0.5700    0.7750    0.3364]) 

    if ~exist(saveDir, 'dir'); mkdir(saveDir); end

    print([saveDir, saveName],'-dpng')
    print([saveDir, saveName],'-deps')
    saveas(gcf, [saveDir, saveName, '.fig'])
end