close all
clear 
figureCount = 1;

load uCTD_SODA_v1.mat
load uCTD_transectStart.mat
load profiles

uCTDtimes = uCTD.time';
uCTDlons = uCTD.lon';
uCTDlats = uCTD.lat';
uCTDtimes = uCTD.time';
startTime = uCTDtimes(uCTD_transectStart(1));
endTime = uCTDtimes(uCTD_transectStart(end));

uCTDlons = uCTDlons(min(uCTD_transectStart):max(uCTD_transectStart));
uCTDlats = uCTDlats(min(uCTD_transectStart):max(uCTD_transectStart));
uCTDtimes = uCTDtimes(min(uCTD_transectStart):max(uCTD_transectStart));

plot_uCTD_met_spatialOverlap = 0;
plot_windMagnitude_check = 0;
plot_shipwinds_SAR = 1;
plot_eraInterim = 0;
mean_eraInterim = 0; 


load met.mat
% for i = 1:56
%     disp(met.README(i))
% end

metTime = met.time;
metLat = met.LA; %LA: NMEA Latitude
metLon = met.LO; %LO: NMEA Longitude
%metTemp = met.TT; %TT: TSG Temperature (deg C)
%metSalt = met.SA; %SA: Salinity (psu)
metAirTemp = met.AT; %AT: Air Temperature (deg C)
%metPres = met.BP; %BP: Barometric Pressure (mb)
%metPAR = met.PA; %PA: Surface PAR (uE/sec/m^2)
metWindSpeed = met.TW; %TW: True Wind Speed (m/s)
metWindDir = met.TI; %TI: True Wind Direction (direction wind is coming from)

met_inTimeRange = zeros(size(metTime));
met_inTimeRange(metTime >= startTime & metTime <= endTime) = 1;
metTime = metTime(met_inTimeRange == 1);
metLat = metLat(met_inTimeRange == 1);
metLon = metLon(met_inTimeRange == 1);
metAirTemp = metAirTemp(met_inTimeRange == 1);

metWindSpeed = metWindSpeed(met_inTimeRange == 1);
metWindDir = metWindDir(met_inTimeRange == 1);

dayFraction_5min = 5 / (60*24); %Use data taken within 5 minutes of the cast

for profNum = 1:length(profiles)% = min(uCTD_transectStart(:)):max(uCTD_transectStart(:))
    profTime = uCTDtimes(profNum);
    met_inRange = find(metTime > (profTime - dayFraction_5min) & metTime < (profTime + dayFraction_5min));
    
    profiles(profNum).metLat = mean(metLat(met_inRange));
    profiles(profNum).metLon = mean(metLon(met_inRange));
    
    profiles(profNum).meanWindSpeed = mean(metWindSpeed(met_inRange));
    profiles(profNum).meanWindDir_compass = mean(metWindDir(met_inRange));
     
    profiles(profNum).meanAirTemp = mean(metAirTemp(met_inRange));
   
end

compassWind = ([profiles.meanWindDir_compass])';
idx = find(compassWind >= 0 & compassWind < 90);
polarCoordWind(idx, 1) = abs(compassWind(idx) - 90);
idx = find(compassWind >=90 & compassWind <= 360);
polarCoordWind(idx,1) = abs(450 - compassWind(idx));

heading = deg2rad(polarCoordWind - 180); 
[u,v] = pol2cart(heading,[profiles.meanWindSpeed]'); %Because wind dir is the direction wind is coming from, and w want the direction wind is blowing toward

for profNum = 1:length(profiles)
    profiles(profNum).u10_met = u(profNum);
    profiles(profNum).v10_met = v(profNum);
end

minlon = min(uCTDlons(:)) - 0.2; 
maxlat = 74.4; 
maxlon = max(uCTDlons(:)) +.4;
minlat = 74.15;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);


if plot_shipwinds_SAR == 1

% Use this code to make the background SAR image. Otherwise just open the image, below 
%     figure(figureCount); figureCount = figureCount + 1;
%     set(gcf, 'Position', [176 117 1074  829])
%     load('MIZ1_201810120308_RS2_SCWA_100m.mat')
%     m_pcolor(LON_grid, LAT_grid, A_grid);
%     shading flat
%     colormap gray
%     caxis([0.05 0.2])
%     hold on
%     freezeColors

    open('oct12_SAR_freezeColors.fig')
    
    load healy_cruisetrack.mat
    shipLons = healy.lon; %shipLons(shipLons < 0) = shipLons + 360;
    shipLats = healy.lat;
    [X_healy, Y_healy] = m_ll2xy(shipLons, shipLats);
    scatter(X_healy, Y_healy, 8, 'r', 'filled')
    hold on
    
    colormap parula
    [uCTDX, uCTDY] = m_ll2xy(uCTDlons, uCTDlats);
    hold on
    scatter(uCTDX, uCTDY, 100, 'w', 'filled')
    scatter(uCTDX, uCTDY, 85, [profiles.meanWindSpeed], 'filled')
    caxis([3, 7.5])
    cb = colorbar;
    ylabel(cb, 'Windspeed [m/sec]', 'FontSize', 12)

    m_quiver(uCTDlons, uCTDlats, u, v, 'color', 'w', 'linewidth', 1); %IMPORTANT to use u', v' so that all inputs are column vectors!

    m_grid
end

%These plots show that the correct met data has been identified using the
%time thresholds
if plot_uCTD_met_spatialOverlap == 1
    figure(figureCount); figureCount = figureCount + 1;
    scatter([profiles.profileNumber],[profiles.lat], 'filled')
    hold on
    scatter([profiles.profileNumber],[profiles.metLat])

    figure(figureCount); figureCount = figureCount + 1;
    scatter([profiles.profileNumber],[profiles.lon], 'filled')
    hold on
    scatter([profiles.profileNumber],[profiles.metLon])
end


%Check that wind speed has been maintained during translation to u, v,
%vectors
if plot_windMagnitude_check == 1
    figure(figureCount); figureCount = figureCount + 1;
    plot(sqrt(u.^2 + v.^2), 'linewidth', 1)
    hold on
    plot([profiles.meanWindSpeed], 'LineStyle', '--')
end

plotEveryLon = 1;
plotEveryLat = 1;
if plot_eraInterim == 1 || mean_eraInterim == 1
    eraTime = double(ncread('eraInterim_oct2018.nc', 'time'))/24 + datenum('1900-01-01 00:00:00');
    eraLat = ncread('eraInterim_oct2018.nc', 'latitude');
    eraLon = ncread('eraInterim_oct2018.nc', 'longitude');
    eraslp = ncread('eraInterim_oct2018.nc', 'msl');
    erau10 = ncread('eraInterim_oct2018.nc', 'u10');
    erav10 = ncread('eraInterim_oct2018.nc', 'v10');

    lwd = ncread('eraInterim_oct2018.nc', 'strd');
    
    [eraLon, eraLat] = meshgrid(eraLon, eraLat);
    eraLon = eraLon';
    eraLat = eraLat';
end

if plot_eraInterim == 1
    oct12_ind = find(eraTime == datenum('October 12 2018 18:00:00'));
    
    open('oct12_SAR_freezeColors.fig')
    hold on

    load healy_cruisetrack.mat
    shipLons = healy.lon; %shipLons(shipLons < 0) = shipLons + 360;
    shipLats = healy.lat;
    [X_healy, Y_healy] = m_ll2xy(shipLons, shipLats);
    scatter(X_healy, Y_healy, 8, 'r', 'filled')
    hold on

    [X_era, Y_era] = m_ll2xy(eraLon, eraLat);
    era_windMagnitude = sqrt(erau10(:, :, oct12_ind) .^2 + erav10(:, :, oct12_ind) .^2);
    scatter(X_era(:), Y_era(:), 100, 'w', 'filled')
    scatter(X_era(:), Y_era(:), 85, era_windMagnitude(:), 'filled')
    caxis([3, 7.5])

    cb = colorbar;
    ylabel(cb, 'Windspeed [m/sec]', 'FontSize', 12)

    m_quiver(eraLon(1:plotEveryLon:end, 1:plotEveryLat:end), eraLat(1:plotEveryLon:end, 1:plotEveryLat:end), erau10(1:plotEveryLon:end, 1:plotEveryLat:end, oct12_ind), erav10(1:plotEveryLon:end, 1:plotEveryLat:end, oct12_ind), 0.2, 'g', 'linewidth', 1)
end

if mean_eraInterim == 1
    eraLat_sub = eraLat(eraLat > minlat & eraLat < maxlat & eraLon > minlon + .2 & eraLon < maxlon - .2);
    eraLon_sub = eraLon(eraLat > minlat & eraLat < maxlat & eraLon > minlon + .2 & eraLon < maxlon - .2);
    [X_erasub, Y_erasub] = m_ll2xy(eraLon_sub, eraLat_sub);
    scatter(X_erasub, Y_erasub, 'r') %Shows which data points are selected

    erau_timeseries = zeros(5, 1);
    erav_timeseries = zeros(5, 1);
    eraLWd_timeseries = zeros(5, 1);
    timestep = [datenum('October 12 2018 00:00:00'), datenum('October 12 2018 00:03:00'), datenum('October 12 2018 06:00:00'), datenum('October 12 2018 12:00:00'), datenum('October 12 2018 18:00:00')];
    for i = 1:length(timestep)
        curInd = find(eraTime == timestep(i));

        cur_erau10_sub = erau10(:, :, curInd); 
        cur_erau10_sub = cur_erau10_sub(eraLat > minlat & eraLat < maxlat & eraLon > minlon + .2 & eraLon < maxlon - .2);
        erau_timeseries(i) = mean(cur_erau10_sub(:));

        cur_erav10_sub = erav10(:, :, curInd);
        cur_erav10_sub = cur_erav10_sub(eraLat > minlat & eraLat < maxlat & eraLon > minlon + .2 & eraLon < maxlon - .2);
        erav_timeseries(i) = mean(cur_erav10_sub(:));
        
        cur_eraLWd_sub = erav10(:, :, curInd);
        cur_eraLWd_sub = cur_eraLWd_sub(eraLat > minlat & eraLat < maxlat & eraLon > minlon + .2 & eraLon < maxlon - .2);
        eraLWd_timeseries(i) = mean(cur_eraLWd_sub(:));


    end
    metTime = met.time;
    metWindSpeed = met.TW; 
    
    plot(metTime(metTime > min(timestep) & metTime < max(timestep)), metWindSpeed(metTime > min(timestep) & metTime < max(timestep)), 'linewidth', 1)
    hold on
    scatter(timestep, sqrt(erau_timeseries .^ 2 + erav_timeseries .^ 2), 'filled')
    lgd = legend('Ship-Measured Wind Speed', 'Mean ERA Interim Wind Speed');
    set(lgd, 'fontsize', 12)
    ylabel('Wind Speed [m/sec]', 'fontsize', 12)
    xlim([min(timestep), max(timestep)])

    datetick('x', 'keeplimits')
    xlabel('Time on October 12 [UTC]', 'fontsize', 12)
    

%         erau10_sub = erau10(eraLat > minlat & eraLat < maxlat & eraLon > minlon + .4 & eraLon < maxlon - .4);
%         erav10_sub = erav10(eraLat > minlat & eraLat < maxlat & eraLon > minlon + .4 & eraLon < maxlon - .4);
%         m_quiver(eraLon(1:plotEveryLon:end, 1:plotEveryLat:end), eraLat(1:plotEveryLon:end, 1:plotEveryLat:end), erau10(1:plotEveryLon:end, 1:plotEveryLat:end, oct12_ind), erav10(1:plotEveryLon:end, 1:plotEveryLat:end, oct12_ind), 0.3, 'r', 'linewidth', 1)
end


