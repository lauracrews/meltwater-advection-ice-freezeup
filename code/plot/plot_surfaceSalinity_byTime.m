% Plot the surface salinity measurements from Healy, Seagliders, uCTD
% divided by time periods

% close all; 
saveFigs = false;

saveDir = [rootPath, 'figures/fig2/'];
saveName = 'SSS_fourTimePeriods_v2';
fntsz = 16; %font size for entire figure. Nice to be able to update here - want big font for paper but smaller is better for presentations

%These data already loaded in run_meltwaterAdvection.m
%Load uCTD and Seaglider profiles
% [profiles, goodProfsMask] = loadProfiles;
% 
% %Load underway data from the Healy system 
% metData = load_correct_underwayData;
% 
% %Load Wave Glider data
% wvdata = loadWaveglider;
% 
% %Load sea ice concentration to add to the map
% load AMSR2_2018.mat

clim_salt = [25.5, 27]; %Color limits for plotting
saltOutcrop = 26; %Salinity outcrop approximating front location
[moorings, ~] = defineSODAconstants;

%Set up map projection
minlon = -150; maxlon = -143; minlat = 72; maxlat = 75.5;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

%Letters to label subplots
alphabet = ('a':'z').';
chars = num2cell(alphabet);
chars = chars.';

figure(1); set(gcf, 'color', 'w', 'pos', [144 99 1371 856])%[144 227 1449 728]) %for 1x4 subplots
for splot = 1:4 %Iterate through the time periods

    switch splot
        case 1            
            startTime = datenum('sept 21 2018');
            endTime = datenum('sept 30 2018')-0.001; %Because this really means sampling until midnight on Sept 29/Sept 30         
            iceDay = datenum('Sept 25 2018'); %Day on which to overlay AMSR2 Ice concentration
            
        case 2
            startTime = datenum('sept 30 2018'); %Start sampling at midnight on Sept 29/Sept 30
            endTime = datenum('oct 6 2018')-0.001;          
            iceDay = datenum('Oct 4 2018');
        case 3
            startTime = datenum('oct 6 2018');
            endTime = datenum('oct 11 2018') -0.001;%+ .2;           
            iceDay = datenum('Oct 8 2018');

        case 4
            %SAR image of ice
%             ax1 = axes;
%             subplot(2, 4, splot + 4, ax1)  
%             load('MIZ1_201810120308_RS2_SCWA_100m.mat')
%             m_pcolor(LON_grid, LAT_grid, A_grid);
%             shading flat
%             colormap(ax1, gray)
%             caxis([0.05 0.2])
%             m_grid('fontsize', fntsz, 'linestyle', 'none')
            startTime = datenum('oct 11 2018');% + .2;
            endTime = datenum('oct 15 2018');            
            iceDay = datenum('Oct 12 2018');
    end   
    
    % Add axes for the ice concentration. 
    ax2 = axes; subplot(2, 4, splot + 4, ax2); hold on
    
    [~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
    curIce = squeeze(AMSR2.SIC(:, :, curIceInd));
    curIce(curIce <= 0) = nan;
    h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
%     set(h2, 'facealpha', .8)
    cmocean('ice')
    m_grid('xticklabels', [], 'yticklabels', [], 'linestyle', 'none') 
%     cb2 = colorbar;
%     cbps = get(cb, 'pos'); cbps(2) = cbps(2) - 0.13;
%     set(cb2, 'fontsize', fntsz, 'location', 'southoutside')
%     set(cb2, 'position', cbps)
%     ylabel(cb2, [' AMSR2 sea ice concentration'], 'fontsize', fntsz)%  
    
    m_gshhs_i('patch',0.6*[1 1 1]); %Draw coastline

    % Add another axis for salinity. New axis is necessary to allow different colormap
    ax1 = axes; subplot(2, 4, splot + 4, ax1)

    %Scatter Healy salinity data to that collected in the time limits
    [~, healyStartInd] = min(abs(metData.times - startTime));
    [~, healyEndInd] = min(abs(metData.times - (endTime)));
    h3 = m_scatter(metData.lons(healyStartInd:healyEndInd), metData.lats(healyStartInd:healyEndInd), 5, metData.SAs(healyStartInd:healyEndInd));
    hold on
    
    % Scatter Wave Glider salinity
    if splot == 1        
        %Calculate absolute salinity
        [SA, ~, ~] = gsw_SA_Sstar_from_SP(wvdata.salts(:, 3), gsw_p_from_z(-wvdata.depths(:,  3), wvdata.lats), wvdata.lons, wvdata.lats);
        h4 = m_scatter(wvdata.lons, wvdata.lats, 60, SA, 'filled');
    end
    
    %Scatter Seaglider profiles used on the heat budget transect
    matName = ['heatbuget_transect', num2str(splot), '_profileNumbers.mat'];
    load([rootPath, 'data/heat_budget_transects/', matName], 'profNums')
    m_scatter(profiles.lons(profNums), profiles.lats(profNums), 120, 'm', 'filled');
    
    %Scatter Seaglider and uCTD salinity in upper five meters    
    profNums = identifyProfileNumbers(profiles, minlon, minlat, maxlon, maxlat, startTime, endTime);    
    
    if splot == 4
        inds = find(strcmp(profiles.dataset(profNums), 'uCTD') & (profiles.lons(profNums) < -146.9));
        m_line([min(profiles.lons(profNums(inds))), max(profiles.lons(profNums(inds)))], ...
            [min(profiles.lats(profNums(inds))), max(profiles.lats(profNums(inds)))], ...
            'color', 'w', 'linewidth', 11)
    
        %Add ship data again because some got covered by the line
        h3 = m_scatter(metData.lons(healyStartInd:healyEndInd), metData.lats(healyStartInd:healyEndInd), 5, metData.SAs(healyStartInd:healyEndInd));

    end
    
    m_scatter(profiles.lons(profNums), profiles.lats(profNums), 60, nanmean(profiles.SA(1:5, profNums)), 'filled');

%     if splot + 4 == 4
%         m_text(-149.5, 75.5, 'SAR Image Taken on Oct 12', 'fontsize', 14, 'fontweight', 'bold')
%         set(ax1, 'pos', ps)
%     end

    %Finish setting up the salinity axis
    caxis(clim_salt) %Constant color axis everywhere
    cmocean('haline')
    m_grid('xtick', [-149:2:-143], 'fontsize', fntsz)
    ps = get(gca, 'pos');
    if splot == 4
        cb = colorbar('Ticks', [25.5:.5:27]);
        set(cb, 'fontsize', fntsz, 'location', 'eastoutside')
        ylabel(cb, ['Absolute salinity (g/kg)'], 'fontsize', fntsz)%
        set(cb, 'pos', [0.9013    0.2426    0.0138    0.2071])
    end
        
       
    if splot == 1 %Find the ship crossing of the meltwater front in the first time period
        startTime_shipFront = datenum('sept 26 2018') + .4;
        endTime_shipFront = datenum('sept 27 2018') - .4;
        [~, startInd] = min(abs(metData.times - startTime_shipFront));
        [~, endInd] = min(abs(metData.times - endTime_shipFront));

        times_met =  metData.times(startInd:endInd); 
        lats_met = metData.lats(startInd:endInd); lons_met = metData.lons(startInd:endInd); 

        %Identify 26 salinity outcrop = front location in ship data
        [~, ind] = min(abs(metData.SAs(startInd:endInd) - saltOutcrop));
        outcrop_ship = [lons_met(ind), lats_met(ind)]; 
        outcropTime_ship = times_met(ind);

        %Add front location from ship to map
        m_scatter(outcrop_ship(1), outcrop_ship(2), 150, '^', 'w', 'filled')
        m_scatter(outcrop_ship(1), outcrop_ship(2), 100, '^', 'k', 'filled')
%         m_text(outcrop_ship(1) + .1, outcrop_ship(2), ['  Healy, ', datestr(outcropTime_ship, 'dd mmm')], 'color', 'k', 'fontsize', fntsz)
        m_text(outcrop_ship(1) + .2, outcrop_ship(2), ['  ', datestr(outcropTime_ship, 'dd mmm')], 'color', 'k', 'fontsize', fntsz)

    elseif splot == 2 || splot == 3 %Find the Seaglider crossings of the front
        
        profNums = find(strcmp(profiles.dataset, 'sg199') & profiles.times >= startTime & profiles.times <= endTime);
        sg_lons = profiles.lons(profNums); sg_lats = profiles.lats(profNums);
        times = profiles.times(profNums);
        salts = nanmean(profiles.SA(1:5, profNums));
        
         %Identify 26 salinity outcrop = front location in Seaglider data
        [~, indFront] = min(abs(salts - saltOutcrop));

        if isempty(indFront) %This should not occur, though the isohaline barely outcrops on Transect 3
            indFront = length(salts); %Defaults to the end of the transect
        end
            
        outcrop_sg = [sg_lons(indFront), sg_lats(indFront)]; 
        outcropTime_sg = times(indFront);
        m_scatter(outcrop_sg(1), outcrop_sg(2), 150, 'w', '^', 'filled')
        m_scatter(outcrop_sg(1), outcrop_sg(2), 100, 'k', '^', 'filled')
            
        if splot == 2
            m_text(outcrop_sg(1) + .2, outcrop_sg(2), ['  ', datestr(outcropTime_sg, 'dd mmm')], 'color', 'k', 'fontsize', fntsz)
        elseif splot == 3 %Plot end of transect as this enttire transect contains meltwater
            m_text(-148.5, 74.25, datestr(outcropTime_sg, 'dd mmm'), 'color', 'k', 'fontsize', fntsz)
        end     
    end
    
    %Make both axes the same size
    set(ax2, 'pos', ps)
    set(ax1, 'pos', ps)

    %Add mooring locations to map
    m_scatter(moorings.all(:, 2), moorings.all(:, 1), 250, 'k', 'p', 'filled')
%     m_scatter(moorings.all(:, 2), moorings.all(:, 1), 200, 'w', 'p', 'filled')

    obsTitle = ['Observations ' datestr(startTime, 'dd mmm'), ' to ', datestr(endTime-1, 'dd mmm')];
    iceTitle = ['AMSR2 sea ice ' datestr(iceDay, 'dd mmm')];
%     title([obsTitle, sprintf(char('\n')), iceTitle], 'fontsize', fntsz)

    text(0.125,0.95,chars{splot + 4},'Units','normalized','FontSize',fntsz+2, 'fontweight', 'bold') %Used splot + 4+4 so that the temperature panels have letter labels a-d

end

if saveFigs 
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end

clearvars -except rootPath ibcao AMSR2 profiles wvdata metData 
