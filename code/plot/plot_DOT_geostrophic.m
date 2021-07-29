%Calculates geostrophic currents from the CryoSat-2 dynamic ocean
%topography (extracted to matlab using extract_DOT.m). Plots DOT along with
%geostrophic current vectors. 
clearvars -except rootPath AMSR2 profiles wvdata metData
close all;

saveFigs = true;
saveDir = [rootPath, 'figures/fig7/'];
saveName = 'DOT_geostrophic';

%Data that had been extracted from individual .asc files by extract_DOT.m
load('allDOT.mat')
defineSODAconstants;

%Geographic range in which to calculate geostrophic currents
minlon = -152; maxlon = -141; minlat = 72; maxlat = 75.7;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

%Time period to use 
startDay = datenum('Sept 4 2018');
endDay = datenum('oct 3 2018');

%Not enough data in this period to be that meaningful
% startDay = datenum('oct 5 2018');
% endDay = datenum('oct 31 2018');

%Identify DOT measurements in the geographic and time ranges, get rid of
%the rest
inRange = zeros(size(times));
inRange(lons >= minlon & lons <= maxlon & lats >= minlat & lats <= maxlat & times >= startDay & times <= endDay) = 1;
lons = lons(inRange == 1);
lats = lats(inRange == 1);
times = times(inRange == 1);
DOT = DOT(inRange == 1);

%Create a regular lat/lon grid, interpolate the DOT data to that grid
[lon_grid, lat_grid] = meshgrid(minlon:0.5:maxlon, minlat:0.15:maxlat);  %At 75 N, 0.15 degrees longitude ~= 16 km, 0.5 degrees longitude ~= 14.5 km
F_DOT = scatteredInterpolant(lons, lats, DOT);
DOT_grid = F_DOT(lon_grid, lat_grid);

%Also interpolate time to that grid, useful for visualizing where the
%satellite was as time progressed
F_times = scatteredInterpolant(lons, lats, times);
times_grid = F_times(lon_grid, lat_grid);

%Begin selecting data that we will use - exlude areas with too large of
%data gaps to fill with linear interpoation
startRow = 4; endRow = 25; startCol = 1; endCol = 22; %Outer bounding box

%Take the subset of data
lon_sub = lon_grid(startRow:endRow, startCol:endCol);
lat_sub = lat_grid(startRow:endRow, startCol:endCol);
dot_sub = DOT_grid(startRow:endRow, startCol:endCol);

%Delineate the area from which we can include some but not all the data due
%to adjacent data holes
midCol = 9; midRow = 3;

%Really, would be better to interpolate this way, but not working
% DOT_sub_new = interp2(LON_sub(~isnan(DOT_sub)), LAT_sub(~isnan(DOT_sub)), DOT_sub(~isnan(DOT_sub)), LON_sub(:), LAT_sub(:));

%Fill missing values by interpolating along lines of constant latitude
sz = size(lon_sub);
for i = 1:sz(1)
    if i <= midRow %We only want to interpolate along part of the row, because the rest of the row doesn't have enough data
        DOTcur = dot_sub(i, midCol:end);  DOTcur = DOTcur';  
        idx = find(~isnan(DOTcur) == 1);
        curlon = lon_sub(i, midCol:end); curlon= curlon';
        dot_sub(i, midCol:end) =  interp1(curlon(idx), DOTcur(idx), curlon); 
        dot_sub(i, 1:midCol) =  nan; 

    else %Interpolate along the whole row
        DOTcur = dot_sub(i, :);  DOTcur = DOTcur';  
        idx = find(~isnan(DOTcur) == 1);
        curlon = lon_sub(i, :); curlon= curlon';
        dot_sub(i, :) =  interp1(curlon(idx), DOTcur(idx), curlon);   
    end
end

%Smooth DOT with a Gaussian filter
dot_sub = imgaussfilt(dot_sub);

%Begin calculating horizontal gradients
sz = size(lon_sub);
deta_dx = nan .* ones(sz); %Change in eta (smoothed DOT) in x
deta_dy = nan .* ones(sz); %Change in eta (smoothed DOT) in y

%Calculate zonal gradients
for i = 1:sz(1)
    if i > midRow
        d_etax = gradient(dot_sub(i, :));    
        dx = m_lldist(lon_sub(i, :), lat_sub(i, :))' .* 1000; %convert km to m
        deta_dx(i, :) = d_etax / dx(1); %dx should be constant at a given latitude
    else %Only want data from some of the current row
        d_etax = gradient(dot_sub(i, midCol+1:end));    
        dx = m_lldist(lon_sub(i, midCol+1:end), lat_sub(i, midCol+1:end))' .* 1000; %convert km to m
        deta_dx(i, midCol+1:end) = d_etax / dx(1); %dx should be constant at a given latitude
    end
end

%Calculate meridional gradients
for i = 1:sz(2)
    if i <= midCol %Only want data from some of the current column
        d_etay = gradient(dot_sub(midRow+1:end, i));
        dy = m_lldist(lon_sub(midRow+1:end, i), lat_sub(midRow+1:end, i)) .* 1000; %convert km to m
        deta_dy(midRow+1:end, i) = d_etay / dy(1); %dy should be constant everywhere
    else 
        d_etay = gradient(dot_sub(:, i));
        dy = m_lldist(lon_sub(:, i), lat_sub(:, i)) .* 1000; %convert km to m
        deta_dy(:, i) = d_etay / dy(1); %dy should be constant everywhere
    end
end
 
%Calculate geostrophic currents
g = 9.81;
ug = -(g./gsw_f(lat_sub)) .* deta_dy;
vg = (g./gsw_f(lat_sub)) .* deta_dx;
speedg = sqrt(ug .^2 + vg .^2); %Current speed

%% Make plots of the raw DOT data and the DOT data interpolated to the new lat/lon grid

figure; set(gcf, 'pos', [712    82   935   866], 'color', 'w')
subplot(2, 2, 1) %Times of raw satellite measurments 
m_scatter(lons, lats, 20, times, 'filled')
caxis([startDay, endDay])
cb = colorbar;
datetick(cb, 'y', 'mmm-dd')

subplot(2, 2, 2) %Raw satellite data
m_scatter(lons, lats, 20, DOT, 'filled')
cb = colorbar;
ylabel(cb, 'Dynamic ocean topography (m)', 'fontsize', 14)
caxis([-0.05, 0.2])

subplot(2, 2, 3) %Times at interpolated positions
m_scatter(lon_grid(:), lat_grid(:), 20, times_grid(:), 'filled')
hold on
m_scatter(lon_sub(~isnan(dot_sub)), lat_sub(~isnan(dot_sub)), 30, 0.3 .* [1 1 1] , 'linewidth', 0.5) 
caxis([startDay, endDay])
cb = colorbar;
datetick(cb, 'y', 'mmm-dd')

subplot(2, 2, 4) %DOT at interpolated positions
m_scatter(lon_grid(:), lat_grid(:), 30, DOT_grid(:), 'filled')
hold on
m_scatter(lon_sub(~isnan(dot_sub)), lat_sub(~isnan(dot_sub)), 30, 0.3 .* [1 1 1] , 'linewidth', 0.5) 
cb = colorbar;
ylabel(cb, 'Dynamic ocean topography (m)', 'fontsize', 14)
caxis([-0.05, 0.2])

%Format plots and add mooring locations
for splot = 1:4
    subplot(2, 2, splot)
    hold on
    m_grid
    m_scatter(moorings(:, 2), moorings(:, 1), 250, 'k', 'p', 'filled')
    m_scatter(moorings(:, 2), moorings(:, 1), 200, 'w', 'p', 'filled')
end

%% Make map with of DOT with ocean currents, ice conditions

%Optionally update map projection to match other remote sensing figures
maxlon = -142; m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

figure; set(gcf, 'pos', [560   263   894   685], 'color', 'w')

%This axis will show the DOT data
ax1 = axes; 
m_contourf(lon_sub, lat_sub, dot_sub, [-0.05:0.01:0.2], 'linestyle', 'none');
hold on

%Add geostrophic veelocity vectors
scaleFac = 1.7;
m_quiver(lon_sub, lat_sub, ug./scaleFac, vg./scaleFac, 0, 'w')%, 'linewidth', 1) %Zero scaling default for m_scatter, but divide velocities by scaleFac so that the arrows are a reasonable length
m_quiver(-145.25, 72.25, .25/scaleFac, 0, 0, 'k', 'linewidth', 1) %Add a vector key to show scale of velocitiess
m_text(-145.75, 72.33, '0.25 m/sec', 'fontsize', 14)

%Point with velocity ~= 0.25, highlight this to compare that vector to the
%0.25 key vector to make sure the key is correct
% speedg(11, 13)
% m_scatter(lon_sub(11, 13), lat_sub(11, 13), 20, 'b')

% Point with velocity ~= 0.1, to check that the vector key is correct
% m_scatter(lon_sub(12, 4), lat_sub(12, 4), 20, 'r')
% speedg(12, 4)
% m_quiver(-145, 72.05, .1/scaleFac, 0, 'k', 'linewidth', 1)

%Finish formatting DOT plot
ps = get(gca, 'pos');
caxis([-0.05, 0.2])
cmocean('speed')
cb = colorbar;
ylabel(cb, ['Dynamic ocean topography (m) collected', sprintf(char('\n')), 'from ', datestr(startDay, 'mmm dd'), ' to ', datestr(endDay, 'mmm dd yyyy')], 'fontsize', 12)
set(cb, 'fontsize', 12, 'position', [0.8,0.55,0.03,0.35])
m_grid('linestyle', 'none', 'fontsize', 12)

%Add AMSR2 ice concentration on a selected day
load AMSR2_2018.mat
iceDay = datenum('Oct 12 2018');    
ax2 = axes; hold on
[~, curIceInd] = min(abs(AMSR2.mattime - iceDay));
curIce = nanmean(AMSR2.SIC(:, :, curIceInd), 3);
curIce(curIce <= 0) = nan;
h2 = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
% set(h2, 'facealpha', .8)
cmocean('ice')
m_grid('Fontsize', 12)%('linestyle', 'none')  
cb2 = colorbar;
ylabel(cb2, [datestr(iceDay, 'mmm dd yyyy'), ' AMSR2 sea ice concentration'], 'fontsize', 12)%  
set(cb2, 'fontsize', 12, 'position', [0.8 0.16 0.03 0.35])

%Make both axes the same size
set(ax2, 'pos', ps)
set(ax1, 'pos', ps)

%Format map, add mooring locations
m_grid( 'fontsize', 12)
m_scatter(moorings(:, 2), moorings(:, 1), 250, 'k', 'p', 'filled')
m_scatter(moorings(:, 2), moorings(:, 1), 200, 'w', 'p', 'filled')
% m_text(moorings(:, 2), moorings(:, 1), {'    SODA-A', '   SODA-B', '   NAV-SW', '   NAV-SE'}, 'fontsize', 14, 'color', 'w')

%Save figures
if saveFigs 
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end
