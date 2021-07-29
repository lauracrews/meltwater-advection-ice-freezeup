%Identifies profiles taken within a specified geographic region and time
%period. Plots those profiles vs. depth as well as the average profile. 

clearvars -except rootPath AMSR2 profiles wvdata metData
close all

saveFigs = true;
saveDir = [rootPath, 'figures/fig4/'];
saveName = 'meltwaterProfiles';

%Heat content integration parameters
integrationDepth = 40;
integrationIso = 1022;
densChange = 0.25;
plotMLD = false; 

integrationType = 2; %Choose if integrating to a constant depth, isopycnal, or mixed layer deptth
switch integrationType
    case 1 %Itegrate to constant depth
        integrationLimit = integrationDepth;
    case 2 %Integrate to constant isopycnal
        integrationLimit = integrationIso;
    case 3 %Integrate to mixed layer depth
        integrationLimit = densChange;
end

%Use data from between these dates
startTime = datenum('sept 22 2018');
endTime = datenum('oct 28 2018');

%Set up map projection for the entire SODA-A to SODA-B area
minlon = -152; maxlon = -142; minlat = 72; maxlat = 75.7; 
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

%Axis limits for plotting
xlimits_temp = [-1.7, 2];
xlimits_salt = [25, 32];
ylimits_profiles = [0, 100];
xlimits_time =  [startTime-1, endTime - 4];
defineSODAconstants;

%Load data structure of all profiles
% profiles = loadProfiles;

%% Optional code to plot sea ice concentration, over which the locations of profiles can be overlaid
figure(1);
set(gcf, 'pos', [560   263   894   685], 'color', 'w')
load AMSR2_2018.mat
iceDay = datenum('Oct 13 2018');    
[~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
curIce = nanmean(AMSR2.SIC(:, :, curIceInd), 3);
curIce(curIce <= 0) = nan;

ax1 = axes;
h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
set(h2, 'facealpha', .8)
cmocean('ice')
caxis([0, 1])
hold on
m_grid
ps = get(gca, 'pos');
cb = colorbar('position', [0.7809    0.5620    0.0301    0.3642]);
ylabel(cb, [' AMSR2 sea ice concentration'], 'fontsize', 12)% 
title(['Sea ice on ', datestr(AMSR2.mattime(curIceInd), 'mmm dd, yyyy')],  'fontsize', 14)% 
set(ax1, 'pos', ps)

%Add mooring locations to map
m_scatter(moorings(:, 2), moorings(:, 1), 250, 'k', 'p', 'filled')
m_scatter(moorings(:, 2), moorings(:, 1), 200, 'w', 'p', 'filled')
%     m_text(moorings(:, 2), moorings(:, 1), {'    SODA-A', '   SODA-B', '   NAV-SW', '   NAV-SE'}, 'fontsize', 14, 'color', 'w')

%% Iterate through combinations of geographic and time bounds. Identify profiles meeting
%those criteria
plotCount = 1;
for group = 8
    switch group
        case 1
            minlon = -148; maxlon = -146; minlat = 73 + 25/60; maxlat = 74.35; %25
            startTime = datenum('sept 22 2018'); endTime = datenum('oct 3 2018');
            col = blue;
        case 2
            minlon = -148; maxlon = -146; minlat = 73 + 25/60; maxlat = 74.35;
            startTime = datenum('oct 5 2018'); endTime = datenum('oct 11 2018');
            col = purple;
        case 3
            minlon = -148; maxlon = -146; minlat = 73 + 25/60; maxlat = 74.35;
            startTime = datenum('oct 11 2018'); endTime = datenum('oct 15 2018'); 
            col = ltblue;
        case 4
            minlon = -147; maxlon = -145.5; minlat = 74.35; maxlat = 74.67;
            startTime = datenum('sept 21 2018'); endTime = datenum('oct 15 2018'); %endTime = datenum('sept 27 2018');
            col = orange;
        case 5
            minlon = -147; maxlon = -145.5; minlat = 74.35; maxlat = 74.67;
            startTime = datenum('sept 27 2018'); endTime = datenum('oct 6 2018');
            col = orange;
        case 6
            minlon = -147; maxlon = -145; minlat = 74.67; maxlat = 75.1;
            startTime = datenum('sept 22 2018'); endTime = datenum('oct 15 2018'); %endTime = datenum('sept 30 2018'); 
            col = red;
         case 7
            minlon = -147; maxlon = -145; minlat = 74.67; maxlat = 75.1;
            startTime = datenum('sept 30 2018'); endTime = datenum('oct 8 2018');
            col = orange;
        case 8 %Remnant ice area
            minlon = -149.5; maxlon = -147; minlat = 72.75; maxlat = 73.25;
            startTime = datenum('sept 26, 2018'); endTime = datenum('sept 27, 2018');
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
        
    %Profiles within the time bounds
    inTimeMask = zeros(size(profiles.times));
    inTimeMask(profiles.times >= startTime & profiles.times < endTime) = 1;

    %Profiles within the geographic bounds
    inRegionMask = zeros(size(profiles.times));    
    inRegionMask(profiles.lats >= minlat & profiles.lats <= maxlat ...
        & profiles.lons >= minlon & profiles.lons <= maxlon) = 1;
    
    %Identify the profiles in this geographic region in this time period
    if group ~=3
        profNums = find(profiles.qualFlag == 1 & inTimeMask == 1 & inRegionMask == 1); 
    else
        vehicleMask = strcmp(profiles.dataset, 'uCTD'); %Want the uCTD profiles only
        profNums = find(goodProfsMask == 1 & inTimeMask == 1 & inRegionMask == 1 & vehicleMask == 1); 
    end
    
    linstl = '-';

    % Plot profiles on single map with sea ice concentration
    figure(1); 
    h(plotCount) = m_scatter(profiles.lons(profNums), profiles.lats(profNums), 55, 'w', 'filled');
    h(plotCount) = m_scatter(profiles.lons(profNums), profiles.lats(profNums), 40, col, 'filled');
    hold on
    
    %Plot this set of profiles on its own map
    figure(2)
    subplot(2, 4, plotCount); set(gcf, 'color', 'w', 'pos', [114 260 1092 688])
    m_scatter(profiles.lons(profNums), profiles.lats(profNums), 15, col);
    hold on
    m_grid
    title([' Profiles from ', datestr(min(profiles.times(profNums)), 'mmm dd'), ' to ',  datestr(max(profiles.times(profNums)), 'mmm dd')])  
    m_scatter(moorings(:, 2), moorings(:, 1), 100, 'k', 'p')   

    %% Plot average temeprature and salinity profiles
   
    if plotMLD %Calculate the mixed layer depth
        densProf = nanmean(profiles.sigthes(:, profNums), 2);
        mldInd = find(densProf - densProf(1) > densChange, 1);
        mld = profiles.z(mldInd);
    end
 
    figure(3); set(gcf, 'color', 'w', 'pos', [560   627   723   321]) 
    subplot(1, 2, 1) %Temperature profile for this region and this time period
    hold on
    for i = 1:length(profNums) 
        plot(profiles.CT(:, profNums(i)), profiles.z, 'LineWidth', 0.5, 'color', [col, 0.15]); %The addition to the color vector makes it transparent
    end 
    
    %Average temperature profile
    plot(nanmean(profiles.CT(:, profNums), 2), profiles.z, 'linestyle', linstl, 'LineWidth', 2, 'color', col);
%     text(0.2, 3*plotCount, ['n = ', num2str(length(curProfiles.lats)), ' profiles'], 'fontsize', 12, 'color', col)%, 'rotation', 90)
    if plotMLD
        line(xlimits_temp, mld .* [1 1], 'color', 'k', 'linestyle', ':'); 
        text(xlimits_temp(2), mld, ['   Mixed', sprintf('\n'), '   layer', sprintf('\n'), '   depth'], 'color', 'k'); 
    end

    subplot(1, 2, 2) %Salinity profile for this region and this time period
    hold on
    for i = 1:length(profNums) 
        plot(profiles.SA(:, profNums(i)), profiles.z, 'LineWidth', 0.5, 'color', [col, 0.15]); %The addition to the color vector makes it transparent        
    end
    
    if plotMLD; line(xlimits_salt, mld .* [1 1], 'color', 'k', 'linestyle', ':'); end

    %Average salinity profile
    plot(nanmean(profiles.SA(:, profNums), 2), profiles.z, 'linestyle', linstl, 'LineWidth', 2, 'color', col)
%     text(28, 3*plotCount, ['n = ', num2str(length(curProfiles.lats)), ' Profiles'], 'fontsize', 12, 'color', col)%, 'rotation', 90)  
            
    plotCount = plotCount+1;
            
end
%%
%Format map plot
figure(1)
m_scatter(moorings(:, 2), moorings(:, 1), 100, 'k', 'p')   
m_text(moorings(:, 2), moorings(:, 1), {'    SODA-A', '   SODA-B', '   NAV-SW', '   NAV-SE'}, 'fontsize', 12)
m_grid

%%

%Format temperature and salinity plots
figure(3)
subplot(1, 2, 1)
freezingPtProf = gsw_CT_freezing(nanmean(profiles.SA(:, profNums), 2), gsw_p_from_z(-profiles.z, 74));
plot(freezingPtProf, profiles.z, 'k', 'linestyle', '--', 'linewidth', 1)
text(-1.6, 25, 'Freezing temperature', 'Rotation', 90, 'color', 'k', 'fontsize', 10)
xlim(xlimits_temp);
ylim(ylimits_profiles)
xlabel(['Conservative temperature (', sprintf(char(176)), 'C)'], 'fontsize', 12, 'FontWeight', 'bold')
ylabel('Depth (m)', 'fontsize', 12, 'FontWeight', 'bold')
set(gca, 'ydir', 'reverse', 'fontsize', 12)
grid on

subplot(1, 2, 2)
xlim(xlimits_salt);
ylim(ylimits_profiles)
grid on
set(gca, 'ydir', 'reverse', 'fontsize', 12)
xlabel('Absolute salinity (g/kg)', 'fontsize', 12, 'FontWeight', 'bold')
ylabel('Depth (m)', 'fontsize', 12, 'FontWeight', 'bold')

%Save figure
if saveFigs   
    close(figure(1)); close(figure(2));
    figure(3)
    
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end

