%Identifies profiles taken within a specified geographic region and time
%period. Plots those profiles vs. depth as well as the average profile. 

% close all

makeCombinedFigure = true;
saveFigs = false;
saveDir = [rootPath, 'figures/fig4/'];
saveName = 'meltwaterProfiles';

fntsz = 14;
densChange = 0.25; %Mixed layer depth definition
plotMLD = true; 

%Set up map projection for the entire SODA-A to SODA-B area
minlon = -152; maxlon = -142; minlat = 72; maxlat = 75.7; 
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

[moorings, ~] = defineSODAconstants;

%% Identify profiles within geographic and time bounds.
minlon = -149.5; maxlon = -147; minlat = 72.75; maxlat = 73.25;
startTime = datenum('sept 26, 2018'); endTime = datenum('sept 27, 2018');
profNums = identifyProfileNumbers(profiles, minlon, minlat, maxlon, maxlat, startTime, endTime);

%% Optional code to plot sea ice concentration, over which the locations of profiles can be overlaid
figure(2);
set(gcf, 'pos', [560   263   894   685], 'color', 'w')
iceDay = datenum('Sep 26 2018');    
[~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
curIce = squeeze(AMSR2.SIC(:, :, curIceInd));
curIce(curIce <= 0) = nan;

ax1 = axes;
h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
cmocean('ice'); caxis([0, 1])
hold on; m_grid
ps = get(gca, 'pos');
cb = colorbar('position', [0.7809    0.5620    0.0301    0.3642]);
ylabel(cb, 'AMSR2 sea ice concentration', 'fontsize', fntsz)
title(['Sea ice on ', datestr(AMSR2.mattime(curIceInd), 'mmm dd, yyyy')],  'fontsize', fntsz)% 
set(ax1, 'pos', ps)

% Plot profiles on single map with sea ice concentration
m_scatter(profiles.lons(profNums), profiles.lats(profNums), 40, 'k', 'filled');
hold on
title([' Profiles from ', datestr(min(profiles.times(profNums)), 'mmm dd'), ' to ',  datestr(max(profiles.times(profNums)), 'mmm dd')])  

%Add mooring locations to map
m_scatter(moorings.all(:, 2), moorings.all(:, 1), 100, 'k', 'p', 'filled')   

%% Plot temeprature and salinity profiles
col = 'k'; linstl = '-';

%Axis limits for plotting
xlimits_temp = [-1.7, 2];
xlimits_salt = [24.85, 32];
ylimits_profiles = [0, 100];

if plotMLD %Calculate the mixed layer depth
    for i = 1:length(profNums)
        densProf = nanmean(profiles.sigthes(:, profNums(i)), 2);
        mldInd = find(densProf - densProf(1) > densChange, 1);
        profiles.mld(profNums(i)) = profiles.z(mldInd);
    end
    mld = nanmean(profiles.mld(profNums));
    [~, mldInd] = min(abs(profiles.z - mld));
end
%%
figure(1); set(gcf, 'color', 'w')%, 'pos', [560   627   723   321]) 
if makeCombinedFigure
    subplot(2, 3, 5)
else
    subplot(1, 2, 1) %Temperature profile for this region and this time period
end
hold on
for i = 1:length(profNums) %plot all temperature profiles for this region and this time period
    l(1) = plot(profiles.CT(:, profNums(i)), profiles.z, 'LineWidth', 0.5, 'color', [col, 0.15]); %The addition to the color vector makes it transparent
end 

%Average temperature profile
l(2) = plot(nanmean(profiles.CT(:, profNums), 2), profiles.z, 'linestyle', linstl, 'LineWidth', 2, 'color', col);
% if plotMLD
%     line(xlimits_temp, mld .* [1 1], 'color', 'k', 'linestyle', ':'); 
%     text(xlimits_temp(2), mld, ['   Mixed', sprintf('\n'), '   layer', sprintf('\n'), '   depth'], 'color', 'k'); 
% end

if plotMLD; l(4) = scatter(nanmean(profiles.CT(mldInd, profNums), 2), mld, 50, 'r', 'filled'); end

if makeCombinedFigure
    subplot(2, 3, 6)
else
    subplot(1, 2, 2)%Salinity profile for this region and this time period
end
hold on
for i = 1:length(profNums) %plot all salinity profiles for this region and this time period
    plot(profiles.SA(:, profNums(i)), profiles.z, 'LineWidth', 0.5, 'color', [col, 0.15]); %The addition to the color vector makes it transparent        
end

%Average salinity profile
plot(nanmean(profiles.SA(:, profNums), 2), profiles.z, 'linestyle', linstl, 'LineWidth', 2, 'color', col)

% if plotMLD; line(xlimits_salt, mld .* [1 1], 'color', 'k', 'linestyle', ':'); end
if plotMLD; scatter(nanmean(profiles.SA(mldInd, profNums), 2), mld, 50, 'r', 'filled'); end

%Format temperature and salinity plots
if makeCombinedFigure
    subplot(2, 3, 5)
else
    subplot(1, 2, 1) %Temperature profile for this region and this time period
end
freezingPtProf = gsw_CT_freezing(nanmean(profiles.SA(:, profNums), 2), gsw_p_from_z(-profiles.z, 74));
l(3) = plot(freezingPtProf, profiles.z, 'k', 'linestyle', '--', 'linewidth', 1);
% text(-1.6, 25, 'Freezing temperature', 'Rotation', 90, 'color', 'k', 'fontsize', 10)
xlim(xlimits_temp);
xlabel(['Conservative temperature (', sprintf(char(176)), 'C)'], 'fontsize', fntsz)%, 'FontWeight', 'bold')
%text(0.2, 3*plotCount, ['n = ', num2str(length(curProfiles.lats)), ' profiles'], 'fontsize', 12, 'color', col)%, 'rotation', 90)
lgd_profs = legend(l(1:4), 'Observed profiles', 'Mean of observations', 'Freezing temperature', 'Mixed layer depth', 'fontsize', fntsz);
set(lgd_profs, 'pos', [0.5434    0.3667    0.1149    0.0830])

if makeCombinedFigure
    subplot(2, 3, 6)
else
    subplot(1, 2, 2)%Salinity profile for this region and this time period
end
    xlim(xlimits_salt);
xlabel('Absolute salinity (g/kg)', 'fontsize', fntsz)%, 'FontWeight', 'bold')
%text(28, 3*plotCount, ['n = ', num2str(length(curProfiles.lats)), ' Profiles'], 'fontsize', 12, 'color', col)%, 'rotation', 90)  


chars = num2cell('a':'f');
for splot = 1:2
    if makeCombinedFigure
        subplot(2, 3, splot + 4)
        text(0.025,0.95,chars{splot+4},'Units','normalized','FontSize',14, 'fontweight', 'bold')
    else
        subplot(1, 2, splot)
    end
    ylim(ylimits_profiles)
    grid on; box on
    set(gca, 'ydir', 'reverse', 'fontsize', fntsz)
    ylabel('Depth (m)', 'fontsize', fntsz)%, 'FontWeight', 'bold')
end

%% Save figure
if saveFigs   
    close(figure(2));
    figure(1)
    
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end

clearvars -except rootPath AMSR2 profiles wvdata metData

