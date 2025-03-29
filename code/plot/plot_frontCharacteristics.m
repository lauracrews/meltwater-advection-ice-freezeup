% Plots crossing of meltwater front by Wave Glider and Seaglider.
% Calculates and plots horizontal buoyancy gradient, vertical buoyancy
% frequency, geostrophic velocity, Richardson number, and Rossby number.
% The subplots for the density and the horizontal buoyancy gradient
% along the Waveglider track are used in the document. Other calculations
% are discussed in the text. 
 
% close all

saveFigs = false;
saveDir = [rootPath, 'figures/fig4/'];
saveName = 'frontCharacteristics';
makeMap = false;

fntsz = 14;
g = 9.81; %Gravity (m/s^2)
rho0 = 1026; %Reference density [kg/m^3]
mld = 15; %Assumed level of no motion for calculating velocity from thermal wind [m] 
[moorings, colors] = defineSODAconstants; %Mooring locations, colors

%Axis limits and colors for plotting
clim_salt = [25.5, 26.5]; clim_temp = [-1.5, 1];
shallowCol = 'k'; deepCol = 0.5 .* [1 1 1]; %Black and gray

% Load and plot data from Wave Glider 153 as this vehicle has salinity at
% 0.2 m as well as 9 m
vehicleInds = find(wvdata.vehicle == 153);
times = wvdata.times(vehicleInds); lons = wvdata.lons(vehicleInds); lats = wvdata.lats(vehicleInds);
temps = wvdata.temps(vehicleInds, 2:3); salts = wvdata.salts(vehicleInds, 2:3); wvdepth = wvdata.depths(vehicleInds, 2:3); %Select the depthsh that have salinity data

%Calculate conservative temperature, absolute salinity, density
[SA, ~, ~] = gsw_SA_Sstar_from_SP(salts, gsw_p_from_z(-wvdepth, repmat(lats, [1, 2])), lons, lats);
CT = gsw_CT_from_t(SA, temps, gsw_p_from_z(-wvdepth,repmat(lats, [1, 2])));

%Select data from Wave Glider crossing the front using time bounds
%These times should match those used in plot_remnantIceSAR.m
startTime = datenum('sept 22 2018') + .07;
endTime = datenum('sept 22 2018') + .55;
[~, startInd] = min(abs(times - startTime));
[~, endInd] = min(abs(times - endTime));

%Data for Wave Glider during front crossing
lons = lons(startInd:endInd); lats = lats(startInd:endInd);
times = times(startInd:endInd);
CT = CT(startInd:endInd, :); SA = SA(startInd:endInd, :);
wvdepth = wvdepth(startInd:endInd, :);

%Convert all to row vectors
SA = SA'; CT = CT'; salts = salts'; temps = temps'; lons = lons'; lats = lats'; wvdepth = wvdepth'; times = times';
rho = gsw_rho(SA, CT, 0);
rho_smooth = movmean(rho, 5, 2); %Smooth the density fields before calculating gradients
%%
delta = m_lldist(lons, lats) .* 1000; %Distance between consecutive measurements. Convert km to m 
dists = cumsum([0; delta])'; %Cumulative displacement along section.
middists = 0.5 .* (dists(1:end-1) + dists(2:end)); %Midpoint coordinates

drho = diff(rho_smooth, 1, 2);
drho_dx = drho ./ repmat(delta', [2, 1]);

%Calculate horizontal buoyancy gradient
M2 = -(g/rho0) * drho_dx;

%% Begin plotting Wave Glider data 

evenTimes = [datenum('22 sept 2018 02:00'):(1/24):datenum('22 sept 2018 13:00')];
dists_atEvenTimes = interp1(times, dists, evenTimes);

subplot(2, 3, 2) %Density along the Wave Glider transect
yyaxis left
plot(dists./1000, rho(1, :), 'color', shallowCol, 'linestyle', ':', 'linewidth', 1)%, 'Marker', 'o', 'MarkerFaceColor', shallowCol)
hold on
plot(dists./1000, rho(2, :), 'color', deepCol, 'linestyle', ':', 'linewidth', 1)%, 'Marker', 'o', 'MarkerFaceColor', deepCol)
h(1) = plot(dists./1000, rho_smooth(1, :), 'color', shallowCol, 'linestyle', '-', 'linewidth', 1, 'Marker', 'o', 'MarkerFaceColor', shallowCol);
h(2) = plot(dists./1000, rho_smooth(2, :), 'color', deepCol, 'linestyle', '-', 'linewidth', 1, 'Marker', 'o', 'MarkerFaceColor', deepCol);
ylabel('Potential Density (kg/m^3)', 'fontsize', fntsz)
ylim([1020.3, 1021.3])
set(gca, 'YColor', 'k');
%%
for i = 1:length(evenTimes)
    text(dists_atEvenTimes(i)./1000, 1021.32, datestr(evenTimes(i), 'hh:MM'), 'fontsize', 10, 'Rotation', 90)
end
title(['Time on 22 Sep', newline, newline], 'FontWeight', 'normal', 'fontsize', fntsz)

yyaxis right %Additional axis to show the density difference between the upper and lower CTDs
h(3) = plot(dists./1000, rho(2, :) - rho(1, :), 'color', colors.red);
set(gca, 'ycolor', colors.red);
% ylabel('\Delta\rho', 'color', colors.red)
ylim([0, 0.3])

lgd = legend(h(1:3), 'At 0.25 m', 'At 9 m', '\Delta\rho');
set(lgd, 'location', 'southeast', 'fontsize', fntsz)

subplot(2, 3, 3) %Horizontal buoyancy gradient
hold on
plot(middists./1000, M2(1, :) ./ 1e-6, 'color', shallowCol, 'linestyle', '-', 'linewidth', 1)
plot(middists./1000, M2(2, :)./ 1e-6 ,'color', deepCol, 'linestyle', '-', 'linewidth', 1)
ylabel('M^2 (1/s^2) x 10^{-6}')
ylim([-1.2, 0.25])

for i = 1:length(evenTimes)
    text(dists_atEvenTimes(i)./1000, .29, datestr(evenTimes(i), 'hh:MM'), 'fontsize', 10, 'Rotation', 90)
end
title(['Time on 22 Sep', newline, newline], 'FontWeight', 'normal', 'fontsize', fntsz)

%%
%Format subplots
chars = num2cell('a':'c');
for splot = 2:3
    subplot(2, 3, splot)
    text(0.025,0.95,chars{splot},'Units','normalized','FontSize',14, 'fontweight', 'bold')
    xlabel('Distance (km)', 'fontsize', fntsz)
    set(gca, 'fontsize', fntsz)
    xlim([0, max(dists./1000)])
    grid on; box on
end


%%
%Save figure
if saveFigs
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end


if makeMap
    figure(2); set(gcf, 'color', 'w', 'pos', [ 560    92   622   856])


    %Set up map projection
    minlon = -148.4; maxlon = -146; minlat = 72.5; maxlat = 73.2;
    m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat+.1]); %Set up map projection

    %Map - even though this is plotted in plot_remnantIceSAR.m, it's useful to
    %plot here too for comparison with the later Seaglider and Healy transect
    %locations

    ax1 = axes; %First axis is for SAR image
%     subplot(4, 3, 1, ax1)
    fileName = 'WG_20180921_RS2_SCWA_100m.mat';
    load(fileName); 
    titleText = ['SAR image from ', datestr([fileName(8:9), ' ', fileName(10:11), ' ', fileName(4:7)], 'mmm dd'),...
        sprintf('\n'), 'Wave Glider data from ', datestr(startTime, 'mmm dd')];
    m_pcolor(LON_grid, LAT_grid, A_grid);
    shading flat
    colormap gray
    caxis([0 0.2])
    hold on
    m_grid('linestyle', 'none')
    ps = get(gca, 'pos');

    ax2 = axes; %This axis is for the surface salinity 
%     subplot(4, 3, 1, ax2)
    % m_line([lons(1), lons(end)], [lats(1), lats(end)], 'color', 'k', 'linewidth', 1, 'linestyle', '--');
    hold on
    h3 = m_scatter(lons, lats, 60, SA(2, :), 'filled');
    caxis(clim_salt)
    cmocean('haline')
    m_grid
    cb = colorbar;
    ylabel(cb, ['S_A (g/kg)'], 'fontsize', fntsz) 
    set(cb, 'location', 'eastoutside');
    title(titleText, 'fontsize', fntsz)

    set(ax1, 'pos', ps)
    set(ax2, 'pos', ps)
end

clearvars -except rootPath AMSR2 profiles wvdata metData