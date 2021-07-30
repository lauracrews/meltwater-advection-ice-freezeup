%Identifies profiles taken within a specified geographic region and time
%period. Plots those profiles vs. depth as well as the average profile. 

clearvars -except rootPath AMSR2 profiles wvdata metData
close all

saveFigs = true;
saveDir = [rootPath, 'figures/fig4/'];
saveName = 'meltwaterProfiles';

densChange = 0.25; %Mixed layer depth definition
plotMLD = false; 

%Set up map projection for the entire SODA-A to SODA-B area
minlon = -152; maxlon = -142; minlat = 72; maxlat = 75.7; 
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

[moorings, ~] = defineSODAconstants;

%% Identify profiles within geographic and time bounds.
minlon = -149.5; maxlon = -147; minlat = 72.75; maxlat = 73.25;
startTime = datenum('sept 26, 2018'); endTime = datenum('sept 27, 2018');
profNums = identifyProfileNumbers(profiles, minlon, minlat, maxlon, maxlat, startTime, endTime);

%% Optional code to plot sea ice concentration, over which the locations of profiles can be overlaid
figure(1);
set(gcf, 'pos', [560   263   894   685], 'color', 'w')
iceDay = datenum('Sep 26 2018');    
[~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
curIce = nanmean(AMSR2.SIC(:, :, curIceInd), 3);
curIce(curIce <= 0) = nan;

ax1 = axes;
h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
cmocean('ice'); caxis([0, 1])
hold on; m_grid
ps = get(gca, 'pos');
cb = colorbar('position', [0.7809    0.5620    0.0301    0.3642]);
ylabel(cb, 'AMSR2 sea ice concentration', 'fontsize', 12)
title(['Sea ice on ', datestr(AMSR2.mattime(curIceInd), 'mmm dd, yyyy')],  'fontsize', 14)% 
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
xlimits_salt = [25, 32];
ylimits_profiles = [0, 100];

if plotMLD %Calculate the mixed layer depth
    densProf = nanmean(profiles.sigthes(:, profNums), 2);
    mldInd = find(densProf - densProf(1) > densChange, 1);
    mld = profiles.z(mldInd);
end

figure(2); set(gcf, 'color', 'w', 'pos', [560   627   723   321]) 
subplot(1, 2, 1) 
hold on
for i = 1:length(profNums) %plot all temperature profiles for this region and this time period
    plot(profiles.CT(:, profNums(i)), profiles.z, 'LineWidth', 0.5, 'color', [col, 0.15]); %The addition to the color vector makes it transparent
end 

%Average temperature profile
plot(nanmean(profiles.CT(:, profNums), 2), profiles.z, 'linestyle', linstl, 'LineWidth', 2, 'color', col);
if plotMLD
    line(xlimits_temp, mld .* [1 1], 'color', 'k', 'linestyle', ':'); 
    text(xlimits_temp(2), mld, ['   Mixed', sprintf('\n'), '   layer', sprintf('\n'), '   depth'], 'color', 'k'); 
end

subplot(1, 2, 2) 
hold on
for i = 1:length(profNums) %plot all salinity profiles for this region and this time period
    plot(profiles.SA(:, profNums(i)), profiles.z, 'LineWidth', 0.5, 'color', [col, 0.15]); %The addition to the color vector makes it transparent        
end

%Average salinity profile
plot(nanmean(profiles.SA(:, profNums), 2), profiles.z, 'linestyle', linstl, 'LineWidth', 2, 'color', col)

if plotMLD; line(xlimits_salt, mld .* [1 1], 'color', 'k', 'linestyle', ':'); end

%Format temperature and salinity plots
subplot(1, 2, 1)
freezingPtProf = gsw_CT_freezing(nanmean(profiles.SA(:, profNums), 2), gsw_p_from_z(-profiles.z, 74));
plot(freezingPtProf, profiles.z, 'k', 'linestyle', '--', 'linewidth', 1)
text(-1.6, 25, 'Freezing temperature', 'Rotation', 90, 'color', 'k', 'fontsize', 10)
xlim(xlimits_temp);
xlabel(['Conservative temperature (', sprintf(char(176)), 'C)'], 'fontsize', 12, 'FontWeight', 'bold')
%text(0.2, 3*plotCount, ['n = ', num2str(length(curProfiles.lats)), ' profiles'], 'fontsize', 12, 'color', col)%, 'rotation', 90)

subplot(1, 2, 2)
xlim(xlimits_salt);
xlabel('Absolute salinity (g/kg)', 'fontsize', 12, 'FontWeight', 'bold')
%text(28, 3*plotCount, ['n = ', num2str(length(curProfiles.lats)), ' Profiles'], 'fontsize', 12, 'color', col)%, 'rotation', 90)  

for splot = 1:2
    subplot(1, 2, splot)
    ylim(ylimits_profiles)
    grid on
    set(gca, 'ydir', 'reverse', 'fontsize', 12)
    ylabel('Depth (m)', 'fontsize', 12, 'FontWeight', 'bold')
end

%% Save figure
if saveFigs   
    close(figure(1));
    figure(2)
    
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end
