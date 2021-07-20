% Plot the surface ocean temperature measurements from Healy, Seagliders, uCTD
% divided by time periods
clearvars -except AMSR2 profiles wvdata metData
close all

saveFigs = true;
saveDir = [userpath, '/meltwaterAdvection/figures/fig2/'];
saveName = 'SST_fourTimePeriods';

%Load uCTD and Seaglider profiles
%profiles = loadProfiles;

%Load underway data from the Healy system 
% metData = load_correct_underwayData;

%Load Wave Glider data
% wvdata = loadWaveglider;

%Load sea ice concentration to add to the map
% load AMSR2_2018.mat

%Set up map projection
minlon = -150; maxlon = -143; minlat = 72; maxlat = 75.5;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

defineSODAconstants;

figure; set(gcf, 'color', 'w', 'pos', [144 227 1449 728]) %for 1x4 subplots
for splot = 1:4 %Iterate through the time periods

    switch splot
        case 1 
            startTime = datenum('sept 21 2018');
            endTime = datenum('sept 30 2018');          
            iceDay = datenum('Sept 25 2018'); %Day on which to overlay AMSR2 Ice concentration
            clim_temp = [-1.5, 1];
            
            %MODIS SST data
            curFile = 'T20182572018264.L3m_8D_NSST_sst_4km.nc'; 
%             curFile = 'T20182652018272.L3m_8D_NSST_sst_4km.nc'; 
        case 2
            startTime = datenum('sept 30 2018');
            endTime = datenum('oct 6 2018');          
            iceDay = datenum('Oct 4 2018');
            clim_temp = [-1.5, 0.5];
            
            %MODIS SST data
            curFile = 'T20182732018280.L3m_8D_NSST_sst_4km.nc';          
        case 3
            startTime = datenum('oct 6 2018');
            endTime = datenum('oct 11 2018');           
            iceDay = datenum('Oct 8 2018');
            clim_temp = [-1.5, 0];
        case 4
            %SAR image of ice - now use AMSR instead for consistency
%             ax1 = axes;
%             subplot(1, 4, splot, ax1)  
%             load('MIZ1_201810120308_RS2_SCWA_100m.mat')
%             m_pcolor(LON_grid, LAT_grid, A_grid);
%             shading flat
%             colormap(ax1, gray)
%             caxis([0.05 0.2])
%             m_grid('fontsize', 12, 'linestyle', 'none')
            startTime = datenum('oct 11 2018');
            endTime = datenum('oct 15 2018');            
            iceDay = datenum('Oct 12 2018');
            clim_temp = [-1.5, -0.25];
    end
    
    if exist([userpath, '/meltwaterAdvection/data/ModisTerra/', curFile], 'file') == 0
        disp(['MODIS-Terra Sea surface temperature data not found. Download the file ', newline,...
            curFile, ' from ', newline, 'https://podaac-tools.jpl.nasa.gov/drive/files/allData/modis/L3/terra/11um/v2019.0/4km/8day/2018', ...
            newline, 'to the directory ~/meltwaterAdvection/data/ModisTerra/'])
        
        return
    end
    
    ax1 = axes; subplot(1, 4, splot, ax1)  
    modisTitle = ''; 
    if splot == 1 || splot == 2 %Plot the MODIS image
        [modisLon, modisLat, sst, ~, ~, modisTitle] = loadMODISsst(curFile, minlon, maxlon, minlat, maxlat);
        h1 = m_pcolor(modisLon, modisLat, sst);
        hold on  
        set(h1, 'facealpha', .9)
    end
   
    %Scatter ship temperatures
    [~, healyStartInd] = min(abs(metData.times - startTime));
    [~, healyEndInd] = min(abs(metData.times - (endTime)));
    m_scatter(metData.lons(healyStartInd:healyEndInd), metData.lats(healyStartInd:healyEndInd), 5, metData.CTs(healyStartInd:healyEndInd));
    hold on
    
    % Scatter Wave Glider temperature
    if splot == 1        
        %Calculate conservative tempereature, absolute salinity
        [SA, ~, ~] = gsw_SA_Sstar_from_SP(wvdata.salts(:, 3), gsw_p_from_z(-wvdata.depths(:,  3), wvdata.lats), wvdata.lons, wvdata.lats);
        CT = gsw_CT_from_t(SA, wvdata.temps(:,  3), gsw_p_from_z(-wvdata.depths(:,  3), wvdata.lats));
        m_scatter(wvdata.lons, wvdata.lats, 60, CT, 'filled');
        
    end
       
    %Scatter Seaglider and uCTD observations in upper five meters
    profNums = identifyProfileNumbers(profiles, minlon, minlat, maxlon, maxlat, startTime, endTime);
    m_scatter(profiles.lons(profNums), profiles.lats(profNums), 60, nanmean(profiles.CT(1:5, profNums)), 'filled');
      
    cmocean('thermal'); caxis(clim_temp); %Set color axis
    ps = get(gca, 'pos');
    cb = colorbar('Ticks', [-1.5, -1, -0.5, 0, 0.5, 1]);
    set(cb, 'fontsize', 12, 'location', 'southoutside')
    ylabel(cb, ['Conservative temperature ', ' (', sprintf(char(176)), 'C)'], 'fontsize', 12)%
    m_grid('xtick', [-149:2:-143], 'fontsize', 12)   

    % Add another axes for the ice concentration. New axis is necessary to allow different colormap
    ax2 = axes; subplot(1, 4, splot, ax2); hold on  
    [~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
    curIce = nanmean(AMSR2.SIC(:, :, curIceInd), 3);
    curIce(curIce <= 0) = nan;
    h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
    set(h2, 'facealpha', .8)
    cmocean('ice')
    m_grid('xticklabels', [], 'yticklabels', [], 'linestyle', 'none')
    cb2 = colorbar;
    cbps = get(cb, 'pos'); cbps(2) = cbps(2) - 0.13;
    set(cb2, 'fontsize', 12, 'location', 'southoutside')
    set(cb2, 'position', cbps)
    ylabel(cb2, [' AMSR2 sea ice concentration'], 'fontsize', 12)%  
    
    %Make both axes the same size
    set(ax2, 'pos', ps)
    set(ax1, 'pos', ps)

    if splot == 1
        %Outline box for zoom in in meltwater area
        minlon_box = -149; maxlon_box = -146.5; minlat_box = 72.5; maxlat_box = 73.2;
        pt1 = [minlon_box, minlat_box]; pt2 = [maxlon_box, minlat_box]; pt3 = [maxlon_box, maxlat_box]; pt4 = [minlon_box, maxlat_box];
        pts = [pt1(1), pt2(1), pt3(1), pt4(1); pt1(2), pt2(2), pt3(2), pt4(2)];
        linpts = [pts'; pts(:, 1)'];
        m_line(linpts(:, 1), linpts(:, 2), 'color', 'c', 'linewidth', 2, 'linestyle', '-')
    end
    
    %Add mooring locations to map
    m_scatter(moorings(:, 2), moorings(:, 1), 250, 'k', 'p', 'filled')
    if splot <= 2; m_scatter(moorings(:, 2), moorings(:, 1), 200, 'w', 'p', 'filled'); end

    m_gshhs_i('patch',0.6*[1 1 1]); %Draw coastline

    obsTitle = ['Observations ' datestr(startTime, 'mmm dd'), ' to ', datestr(endTime, 'mmm dd')];
    iceTitle = ['AMSR2 sea ice ' datestr(iceDay, 'mmm dd')];
    
    title([modisTitle, sprintf('\n'), obsTitle, sprintf('\n'), iceTitle], 'fontsize',14, 'fontweight', 'bold') 
end

if saveFigs 
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end