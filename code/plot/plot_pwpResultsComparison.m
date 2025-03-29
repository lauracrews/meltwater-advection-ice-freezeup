%Plots initial temperature and salinity profiles as well as the 
%resulting PWP profiles at time of freezing for specified regions and time
%periods. Uses a specified example profile fpr each region/time

close all; 
clearvars -except rootPath AMSR2 profiles wvdata metData saveFigs

saveFigs = true;

saveDir = [rootPath, 'figures/fig10/'];
saveName = 'northPWPprofiles';

plotMLD = false;
plotIntegrationDepth = false;

[moorings, colors] = defineSODAconstants;

%Set up map projection for the entire SODA-A to SODA-B area
minlat = 72.95; maxlat = 75.5; minlon = -150; maxlon =  -144;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

%Axis limits for plotting
xlimits_temp = [-1.7, .5];
xlimits_salt = [25.5, 30];
ylimits_profiles = [0, 50];

figure; set(gcf, 'pos', [31 169 1170 690], 'color', 'w');
subplot(1, 3, 1); hold on %Map
m_grid

%Add mooring locations to map
m_scatter(moorings.all(:, 2), moorings.all(:, 1), 250, 'k', 'p', 'filled')

%%
for group = [2, 1]
    %%
    iceDay = [];
    switch group       
         case 1
            col = colors.red;
            profNum = 531;
        case 2
            col = .5 .* [1 1  1];
            profNum = 682; 
            iceDay = datenum('Oct 12 2018');
    end
          
%%   Plot temperature profile 
    subplot(1, 3, 2); hold on
    a(1) = plot(profiles.CT(:, profNum), profiles.z, 'linewidth', 2, 'color', col, 'linestyle', ':');
    a(2) = plot(profiles.pwpFreezeCTprof(:, profNum), profiles.pwpz, 'linewidth', 1, 'color', col);
    freezingPtProf = gsw_CT_freezing(profiles.pwpFreezeSAprof(:, profNum), gsw_p_from_z(-profiles.pwpz, profiles.lats(profNum)));
    
    if group == 2
    fpp = plot(freezingPtProf, profiles.pwpz, 'linewidth', 1, 'linestyle', '--', 'color', 'k'); a(3) = fpp;
    end
%     fpt = text(-1.6, 14, 'Freezing temperature', 'Rotation', 90, 'color', 'k', 'fontsize', 12);
    xlim(xlimits_temp); 
    xlabel(['Conservative temperature (', sprintf(char(176)), 'C)'], 'fontsize', 12, 'FontWeight', 'bold')
 
%%   Plot salinity profile     
    subplot(1, 3, 3); hold on
    a1 = plot(profiles.SA(:, profNum), profiles.z, 'linewidth', 2, 'color', col, 'linestyle', ':');
    a2b = plot(profiles.pwpFreezeSAprof(:, profNum), profiles.pwpz, 'linewidth', 1, 'color', col);
    xlim(xlimits_salt);
    xlabel(['Absolute salinity (g/kg)'], 'fontsize', 12, 'FontWeight', 'bold')

    %% Find the depth of various isopycnals and add horizonal lines at those depths
    if plotIntegrationDepth 
        subplot(1, 3, 2)
        d(2) = line(xlimits_temp, profiles.integrationDepths(profNum) .* [1 1],  'color', 'k', 'linestyle', ':');
        d(5) = text(xlimits_temp(2), profiles.integrationDepths(profNum), ['Integration', sprintf('\n'), 'depth'], 'color', 'k');

        subplot(1, 3, 3)
        d(8) = line(xlimits_salt, profiles.integrationDepths(profNum) .* [1 1],  'color', 'k', 'linestyle', ':');
    end   
    
    if  plotMLD & ~isnan(profiles.pwpFreezeTime(profNum)) 
        densChange = 0.01;
        densProf0 = profiles.sigthes(:, profNum);
        depthInd0 = find(densProf0 - densProf0(1) > densChange, 1);

        densProff = gsw_rho(profiles.pwpFreezeSAprof(:, profNum), profiles.pwpFreezeCTprof(:, profNum), 0);
        depthIndf = find(densProff - densProff(1) > densChange, 1);

        subplot(1, 3, 2)
        b(1) = line(xlimits_temp, profiles.z(depthInd0) .* [1 1], 'color', 'k', 'linestyle', ':'); 
        b(2) = text(xlimits_temp(2), profiles.z(depthInd0), ['  \Delta \rho_0 = ', num2str(densChange)], 'color', col); 

        b(3) = line(xlimits_temp, profiles.pwpz(depthIndf) .* [1 1], 'color', 'k', 'linestyle', ':'); 
        b(4) = text(xlimits_temp(2), profiles.pwpz(depthIndf), ['  \Delta \rho_f = ', num2str(densChange)], 'color', 'k'); 

        subplot(1, 3, 3)
        b(5) = line(xlimits_salt, profiles.z(depthInd0) .* [1 1], 'color', 'k', 'linestyle', ':'); 
        b(6) = line(xlimits_salt, profiles.pwpz(depthIndf) .* [1 1], 'color', 'k', 'linestyle', ':'); 

        %%%%
        densChange = 0.1;
        densProf0 = profiles.sigthes(:, profNum);
        depthInd0 = find(densProf0 - densProf0(1) > densChange, 1);

        densProff = gsw_rho(profiles.pwpFreezeSAprof(:, profNum), profiles.pwpFreezeCTprof(:, profNum), 0);
        depthIndf = find(densProff - densProff(1) > densChange, 1);

        subplot(1, 3, 2)
        b(7) = line(xlimits_temp, profiles.z(depthInd0) .* [1 1], 'color', 'k', 'linestyle', ':'); 
        b(8) = text(xlimits_temp(2), profiles.z(depthInd0), ['  \Delta \rho_0 = ', num2str(densChange)], 'color', col); 

        b(9) = line(xlimits_temp, profiles.pwpz(depthIndf) .* [1 1], 'color', 'k', 'linestyle', ':'); 
        b(10) = text(xlimits_temp(2), profiles.pwpz(depthIndf), ['  \Delta \rho_f = ', num2str(densChange)], 'color', 'k'); 

        subplot(1, 3, 3)
        b(11) = line(xlimits_salt, profiles.z(depthInd0) .* [1 1], 'color', 'k', 'linestyle', ':'); 
        b(12) = line(xlimits_salt, profiles.pwpz(depthIndf) .* [1 1], 'color', 'k', 'linestyle', ':'); 
    end
    
    %Make legend now that we are finished adding to the T-S plot
    
    lgda = legend(a(1:3), 'Initial observed profile', 'PWP profile at freeze up', 'Freezing temperature');
    set(lgda, 'fontsize', 12, 'location', 'south');
    
    % Add profile location to map
    subplot(1, 3, 1) %(2, 5, 5)
    m_scatter(profiles.lons(profNum), profiles.lats(profNum),  60, col, 'filled'); %Dot to mark the current location
%         g(2) = m_scatter(profiles.lons(profNum), profiles.lats(profNum),  30, 'w', 'filled');  
    m_text(profiles.lons(profNum)-3.75, profiles.lats(profNum)+.08, ['Initial profile: ', datestr(profiles.times(profNum), 'dd mmm hh:MM')], 'color', col, 'fontsize', 11, 'FontWeight', 'bold')
    m_text(profiles.lons(profNum)-3.75, profiles.lats(profNum)-.03, ['Freeze up: ', datestr(profiles.pwpFreezeTime(profNum), 'dd mmm hh:MM')], 'color', col, 'fontsize', 11, 'FontWeight', 'bold')
    
    %Add ice to the map, if a day for ice wasa defined at the start of the loop
    if ~isempty(iceDay)
        load AMSR2_2018.mat
        [~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
        curIce = nanmean(AMSR2.SIC(:, :, curIceInd), 3); curIce(curIce <= 0) = nan;
        m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
        cmocean('ice')
        cb = colorbar; set(cb, 'fontsize', 12, 'location', 'southoutside')
        ylabel(cb, ['AMSR2 sea ice concentration', sprintf('\n'), ' on ', datestr(iceDay, 'dd mmm')], 'fontsize', 12)%  
    end

end
  %%
%Format the temperature and salinity profiles
for splot = 2:3
    subplot(1, 3, splot)
    grid on; box on
    ylim(ylimits_profiles)
    if splot == 2; ylabel('Depth (m)', 'fontsize', 12, 'FontWeight', 'bold'); end
    set(gca, 'ydir', 'reverse', 'fontsize', 12)
end
    
%Save figure
if saveFigs    
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
    print([saveDir, saveName],'-deps')

end

