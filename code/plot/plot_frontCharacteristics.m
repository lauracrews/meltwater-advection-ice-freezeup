% Plots crossing of meltwater front by Wave Glider and Seaglider.
% Calculates and plots horizontal buoyancy gradient, vertical buoyancy
% frequency, geostrophic velocity, Richardson number, and Rossby number.
% The subplots for the density and the horizontal buoyancy gradient
% along the Waveglider track are used in the document. Other calculations
% are discussed in the text. 
 
close all
clearvars -except rootPath AMSR2 profiles wvdata metData

saveFigs = true;
saveDir = [rootPath, 'figures/fig4/'];
saveName = 'frontCharacteristics';

plotAdditionalCrossings = false; %Switch to plot additional front crossings by Healy and Seaglider (not included in paper)

g = 9.81; %Gravity [m/sec^2]
rho0 = 1026; %Reference density [kg/m^3]
mld = 15; %Assumed level of no motion for calculating velocity from thermal wind [m] 
[moorings, colors] = defineSODAconstants; %Mooring locations, colors

%Axis limits and colors for plotting
clim_salt = [25.5, 26.5]; clim_temp = [-1.5, 1];
shallowCol = 'k'; deepCol = 0.5 .* [1 1 1]; %Black and gray

% Load and plot data from Wave Glider 153 as this vehicle has salinity at
% 0.2 m as well as 9 m
wvdata = loadWaveglider;
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
CT = CT(startInd:endInd, :); SA = SA(startInd:endInd, :);
wvdepth = wvdepth(startInd:endInd, :);

%Convert all to row vectors
SA = SA'; CT = CT'; salts = salts'; temps = temps'; lons = lons'; lats = lats'; wvdepth = wvdepth';
rho = gsw_rho(SA, CT, 0);
rho_smooth = movmean(rho, 5, 2); %Smooth the density fields before calculating gradients

%Calculate average density and depth between points to approximate values at midpoints 
rho_mid = .5 .* (rho(:, 1:end-1) +rho(:, 2:end));
wvdepth_mid = .5 .* (wvdepth(:, 1:end-1) +wvdepth(:, 2:end));

%Calculate buoyancy frequency
N2 = (g./rho0) .* (rho_mid(2,:) - rho_mid(1,:)) ./ (wvdepth_mid(2,:) - wvdepth_mid(1,:));

delta = m_lldist(lons, lats) .* 1000; %Distance between consecutive measurements. Convert km to m 
dists = cumsum([0; delta])'; %Cumulative displacement along section.
middists = 0.5 .* (dists(1:end-1) + dists(2:end)); %Midpoint coordinates

%Calculate rho derivatives in x and y directions from the along-track data
midlats = nan .* ones(size(lats(1:end-1)));
midlons = nan .* ones(size(lats(1:end-1)));
drho_dx = nan .* ones(size(rho(:, 1:end-1)));
drho_dy = nan .* ones(size(rho(:, 1:end-1)));
dx = nan .* ones(size(lats(1:end-1)));
dy = nan .* ones(size(lats(1:end-1)));
for i = 1:length(lats) - 1
    midlats(i) = mean(lats(1:i+1));
    midlons(i) = mean(lons(1:i+1));
    
    [delta, a12, ~] = m_idist(lons(i), lats(i), lons(i+1), lats(i+1)); %Distance and bearing between consecutive points
    theta = mod(450 - a12, 360); %Convert compass coordinate to cartesian
    drho = rho_smooth(:, i+1) - rho_smooth(:, i); %Change in rho between consecutive points
    
    dy(i) = delta * sind(theta); %y-dir displacement
    dx(i) = delta * cosd(theta); %x-dir displacement
    drho_dy(:, i) = (drho/delta) * sind(theta); %y-dir gradient in rho
    drho_dx(:, i) = (drho/delta) * cosd(theta); %x-dir gradient in rho
end

%Coriolis at midpointsf
midf = gsw_f(midlats); 

%Calculate horizontal buoyancy gradients
M2_x = drho_dx * g / rho0;
M2_y = drho_dy * g / rho0;
M2 = sqrt(M2_x.^2 + M2_y.^2);

%Calculate Richardon number using upper and lower horizontal buoyancy gradients
Ri(1, :) = midf.^2 .* N2 ./ (M2(1, :).^2);
Ri(2, :) = midf.^2 .* N2 ./ (M2(2, :).^2);

%Calculate vertical velocity shear using geostrophic balance
dv_dz = -g./(rho0 .* repmat(midf, [2,1])) .* drho_dx;
du_dz = g./(rho0 .* repmat(midf, [2,1])) .* drho_dy;

%Calculate geostrophic velocities assuming constant shear, no motion at
%level of mld
ug = du_dz .* (mld-wvdepth_mid);
vg = dv_dz .* (mld-wvdepth_mid);
speed = sqrt(ug.^2 + vg.^2);

% Calculate horizontal velocity gradients, relative vorticity, Rossby number
[dv_dx, ~] = gradient(vg, repmat(cumsum([0, dx]), 2, 1), 1);
[du_dy, ~] = gradient(ug, repmat(cumsum([0, dy]), 2, 1), 1);
Ro = (dv_dx - du_dy) ./ repmat(midf, 2, 1);

%% Begin plotting Wave Glider data 
figure(1); set(gcf, 'color', 'w', 'pos', [ 560    92   622   856])


%Set up map projection
minlon = -148.4; maxlon = -146; minlat = 72.5; maxlat = 73.2;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat+.1]); %Set up map projection

%Map - even though this is plotted in plot_remnantIceSAR.m, it's useful to
%plot here too for comparison with the later Seaglider and Healy transect
%locations

ax1 = axes; %First axis is for SAR image
subplot(4, 3, 1, ax1)
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
subplot(4, 3, 1, ax2)
% m_line([lons(1), lons(end)], [lats(1), lats(end)], 'color', 'k', 'linewidth', 1, 'linestyle', '--');
hold on
h3 = m_scatter(lons, lats, 60, SA(2, :), 'filled');
caxis(clim_salt)
cmocean('haline')
m_grid
cb = colorbar;
ylabel(cb, ['S_A [g/kg]'], 'fontsize', 14) 
set(cb, 'location', 'eastoutside');
title(titleText, 'fontsize', 12)

set(ax1, 'pos', ps)
set(ax2, 'pos', ps)

subplot(4, 2, 3) %Density along the Wave Glider transect
yyaxis left
plot(dists./1000, rho(1, :), 'color', shallowCol, 'linestyle', ':', 'linewidth', 1)%, 'Marker', 'o', 'MarkerFaceColor', shallowCol)
hold on
plot(dists./1000, rho(2, :), 'color', deepCol, 'linestyle', ':', 'linewidth', 1)%, 'Marker', 'o', 'MarkerFaceColor', deepCol)
h(1) = plot(dists./1000, rho_smooth(1, :), 'color', shallowCol, 'linestyle', '-', 'linewidth', 1, 'Marker', 'o', 'MarkerFaceColor', shallowCol);
h(2) = plot(dists./1000, rho_smooth(2, :), 'color', deepCol, 'linestyle', '-', 'linewidth', 1, 'Marker', 'o', 'MarkerFaceColor', deepCol);
ylabel('Potential Density')
ylim([1020.3, 1021.3])
set(gca, 'YColor', 'k');

yyaxis right %Additional axis to show the density difference between the upper and lower CTDs
h(3) = plot(dists./1000, rho(2, :) - rho(1, :), 'color', colors.red);
set(gca, 'ycolor', colors.red);
% ylabel('\Delta\rho', 'color', colors.red)
ylim([0, 0.3])

lgd = legend(h(1:3), 'At 0.25 m', 'At 9 m', '\Delta\rho');
set(lgd, 'location', 'southeast', 'fontsize', 12)

subplot(4, 2, 4) %Meridional geostrophic velocity
hold on
% plot(middists./1000, ug(1,:));
% plot(middists./1000, ug(2,:));
% plot(middists./1000, sqrt(ug(1, :).^2 + vg(1, :).^2), 'color', shallowCol, 'linestyle', '-', 'linewidth', 1)
% plot(middists./1000, sqrt(ug(2, :).^2 + vg(2, :).^2), 'color', deepCol, 'linestyle', '-', 'linewidth', 1);
plot(middists./1000, vg(1, :), 'color', shallowCol, 'linestyle', '-', 'linewidth', 1)
plot(middists./1000, vg(2, :), 'color', deepCol, 'linestyle', '-', 'linewidth', 1);
ylabel('v_{geostrophic} [m/sec]')

subplot(4, 2, 5) %Horizontal buoyancy gradient
hold on
plot(middists./1000, M2(1, :), 'color', shallowCol, 'linestyle', '-', 'linewidth', 1)
plot(middists./1000, M2(2, :),'color', deepCol, 'linestyle', '-', 'linewidth', 1)
ylabel('M^2 [sec^{-2}]')

subplot(4, 2, 6) %Rossby number
plot(middists./1000, Ro(1, :), 'color', shallowCol, 'linestyle', '-', 'linewidth', 1)
hold on
plot(middists./1000, Ro(2, :), 'color', deepCol, 'linestyle', '-', 'linewidth', 1)
ylabel('M^2 [sec^{-2}]')
ylabel('Ro [\zeta / f]', 'fontsize', 14)

subplot(4, 2, 7) %Buoyancy frequency
plot(middists./1000, N2, 'color', shallowCol, 'linestyle', '-', 'linewidth', 1)
hold on
ylabel('N^2 [sec^{-2}]')

subplot(4, 2, 8) %Richardson number
plot(middists./1000, Ri(1, :), 'color', shallowCol, 'linestyle', '-', 'linewidth', 1)
hold on
plot(middists./1000, Ri(2, :), 'color', deepCol, 'linestyle', '-', 'linewidth', 1)
ylabel('Ri [f^2 N^2 M^{-4}]')
ylim([0, 20])

%Format subplots
for splot = 3:8
    subplot(4, 2, splot)
    xlabel('Distance [km]')
    xlim([0, max(dists./1000)])
    grid on
end

%Save figure
if saveFigs
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end

%Decide to continue on and plot later Healy and Seaglider data
if ~plotAdditionalCrossings
    return
end

%% Density gradient from Healy underway system
metData = load_correct_underwayData;
rho = gsw_rho(metData.SAs, metData.CTs, 0);

%Time range of ship data in which the front crossing occurred 
startTime_shipFront = datenum('sept 26 2018') + .46;
endTime_shipFront = datenum('sept 26 2018') + .54;

%Subset data to the front crossing
[~, startInd]= min(abs(metData.times - startTime_shipFront));
[~, endInd] = min(abs(metData.times - endTime_shipFront));
lats = metData.lats(startInd:endInd); lons = metData.lons(startInd:endInd);
salts = metData.SAs(startInd:endInd); rho= rho(startInd:endInd);
rho_mid = .5 .* (rho(:, 1:end-1) +rho(:, 2:end)); %Approximate density at midpoints by averaging

delta = m_lldist(lons', lats') .* 1000; %Distance betweeen measurements along section. Convert km to m 
dists = cumsum([0; delta])'; %Displacement along section
middists = 0.5 .* (dists(1:end-1) + dists(2:end)); %Midpoint coordinates
% figure; hist(delta); %To  see horizontal resolution of Healy samples

rho_smooth = movmean(rho, 3, 2); %Smooth density data before calculating gradients

%Calculate rho derivatives in x and y directions from the along-track data
midlats = nan .* ones(size(lats(1:end-1)));
midlons = nan .* ones(size(lats(1:end-1)));
drho_dx = nan .* ones(size(rho(1:end-1)));
drho_dy = nan .* ones(size(rho(1:end-1)));
dx = nan .* ones(size(lats(1:end-1)));
dy = nan .* ones(size(lats(1:end-1)));
for i = 1:length(lats) - 1
    midlats(i) = mean(lats(1:i+1));
    midlons(i) = mean(lons(1:i+1));
    
    [delta, a12, ~] = m_idist(lons(i), lats(i), lons(i+1), lats(i+1));
    theta = mod(450 - a12, 360); %Convert compass coordinate to cartesian
    drho = rho_smooth(i+1) - rho_smooth(i);
    drho_dy(i) = (drho/delta) * sind(theta);
    drho_dx(i) = (drho/delta) * cosd(theta);
    dy(i) = delta * sind(theta);
    dx(i) = delta * cosd(theta);
end

M2_x = drho_dx * g / rho0;
M2_y = drho_dy * g / rho0;
M2 = sqrt(M2_x.^2 + M2_y.^2);

subplot(4, 3, 2) %Map with front crossing
minlon = -148; maxlon = -146; minlat = 73.1; maxlat = 74;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);
hold on
h3 = m_scatter(lons, lats, 60, salts, 'filled');
caxis(clim_salt)
cmocean('haline')
m_scatter(moorings.all(2:3, 2), moorings.all(2:3, 1), 250, 'k', 'p', 'filled')
m_scatter(moorings.all(2:3, 2), moorings.all(2:3, 1), 200, 'w', 'p', 'filled')
m_grid
% cb = colorbar;
% ylabel(cb, ['SA [g/kg]'], 'fontsize', 14) 

subplot(4, 2, 3) %Density along section
plot(dists./1000, rho, 'color', shallowCol, 'linestyle', '-.', 'linewidth', 1)%, 'Marker', 'o', 'MarkerFaceColor', shallowCol)
hold on
ylabel('Potential Density')
% legend('Surface', 'At 9 m')
% ylim([1021.3, 1021.7])

subplot(4, 2, 5) %Horizontal buoyancy gradient
hold on
plot(middists./1000, M2, 'color', shallowCol, 'linestyle', '-.', 'linewidth', 1)
ylabel('M^2 [sec^{-2}]')

%% Data and plots from Seaglider section

z = [1, 9]; %Depths to plot
mldThreshold = 0.25; %Density change threshold relative to surface for plotting mixed layer depth

% Load data structure of all profiles
[profiles, goodProfsMask] = loadProfiles;

startTime = datenum('sept 30 2018 06:00:00')-.5;
endTime = datenum('oct 5 2018 20:00:00');   
inTimeMask = zeros(size(profiles.times));
    inTimeMask(profiles.times >= startTime & profiles.times < endTime) = 1;

%In the geographic bounds
minlon = -147.75; maxlon = -146.25; minlat = 73.5; maxlat = 74.1;
inRegionMask = zeros(size(profiles.times));    
inRegionMask(profiles.lats >= minlat & profiles.lats <= maxlat ...
    & profiles.lons >= minlon & profiles.lons <= maxlon) = 1;

vehicleMask = strcmp(profiles.dataset, 'sg199');

profNums = find(goodProfsMask == 1 & inTimeMask == 1 & inRegionMask == 1 & vehicleMask == 1); 
    
%Calculate mixed layer depth
for i = 1:length(profiles.lats)
    mldind = find(profiles.sigthes(:, i) - profiles.sigthes(1, i) >= mldThreshold, 1, 'first');
    profiles.mld(i) = profiles.z(mldind);
end

temps = profiles.CT(:, profNums); salts = profiles.SA(:, profNums); sigthes = profiles.sigthes(:, profNums); mlds = profiles.mld(profNums);
lats = profiles.lats(profNums); lons = profiles.lons(profNums);

% %Update arrays to be sorted by latitude
[~, idx] = sort(lats);
temps = temps(:, idx); salts = salts(:, idx); sigthes = sigthes(:, idx); 
lats = lats(idx); lons = lons(idx); mlds = mlds(idx); 

delta = m_lldist(lons', lats') .* 1000; %Distance between consecutive measurements. Convert km to m 
dists = cumsum([0; delta])'; %Displacement along section.
middists = 0.5 .* (dists(1:end-1) + dists(2:end));

rho = sigthes(1:30, :);
rho_smooth = movmean(rho, 5, 2);
rho_mid = .5 .* (rho(:, 1:end-1) + rho(:, 2:end));
N2 = (g./rho0) .* (rho_mid(z(2),:) - rho_mid(z(1),:)) ./ (z(2) - z(1));

midlats = nan .* ones(size(lats(1:end-1)));
midlons = nan .* ones(size(lats(1:end-1)));
drho_dx = nan .* ones(size(rho(:, 1:end-1)));
drho_dy = nan .* ones(size(rho(:, 1:end-1)));
dx = nan .* ones(size(lats(1:end-1)));
dy = nan .* ones(size(lats(1:end-1)));
for i = 1:length(lats) - 1
    midlats(i) = mean(lats(1:i+1));
    midlons(i) = mean(lons(1:i+1));
    
    [delta, a12, ~] = m_idist(lons(i), lats(i), lons(i+1), lats(i+1));
    theta = mod(450 - a12, 360); %Convert compass coordinate to cartesian
    drho = rho_smooth(:, i+1) - rho_smooth(:, i);
    drho_dy(:, i) = (drho/delta) * sind(theta);
    drho_dx(:, i) = (drho/delta) * cosd(theta);
    dy(i) = delta * sind(theta);
    dx(i) = delta * cosd(theta);
end

midf = gsw_f(midlats); 

%Calculate horizontal buoyancy gradients
M2_x = drho_dx * g / rho0;
M2_y = drho_dy * g / rho0;
M2 = sqrt(M2_x.^2 + M2_y.^2);

clear Ri
Ri(1, :) = midf.^2 .* N2 ./ (M2(z(1), :).^2);
Ri(2, :) = midf.^2 .* N2 ./ (M2(z(2), :).^2);

dv_dz = -1./repmat(midf, [size(M2_x, 1), 1]) .* M2_x;
du_dz = 1./repmat(midf, [size(M2_x, 1), 1]) .* M2_y;

ug = nan.* ones(size(M2));
vg = nan.* ones(size(M2));
intDepthInd = 30;
for k = intDepthInd:-1:1
    for i = 1:length(midlats)
        ug(k, i) = trapz(du_dz(k:intDepthInd, i));
        vg(k, i) = trapz(dv_dz(k:intDepthInd, i));
    end
end
[du_dy, ~] = gradient(ug, cumsum([0, dy])', profiles.z);
[dv_dx, ~] = gradient(vg, cumsum([0, dx])', profiles.z);
Ro = (dv_dx - du_dy) ./ repmat(midf, [size(M2_x, 1), 1]);

subplot(4, 3, 3)
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat+.1]); %Set up map projection
hold on
h3 = m_scatter(lons, lats, 60, salts(8, :), 'filled');
caxis(clim_salt)
cmocean('haline')
m_grid
title(['Seaglider data from ', datestr(startTime, 'mmm dd'), ' to ', datestr(endTime, 'mmm dd')], 'fontsize', 12)
% cb = colorbar;
% ylabel(cb, ['SA [g/kg]'], 'fontsize', 14) 

subplot(4, 2, 3) %Density along section
plot(dists./1000, rho(z(1), :), 'color', shallowCol, 'linestyle', ':', 'linewidth', 1)
hold on
plot(dists./1000, rho(z(2), :), 'color', deepCol, 'linestyle', ':', 'linewidth', 1)
h(1) = plot(dists./1000, rho_smooth(z(1), :), 'color', shallowCol, 'linestyle', '-', 'linewidth', 1, 'Marker', 'o', 'MarkerFaceColor', shallowCol);
h(2) = plot(dists./1000, rho_smooth(z(2), :), 'color', deepCol, 'linestyle', '-', 'linewidth', 1, 'Marker', 'o', 'MarkerFaceColor', deepCol);

ylabel('Potential Density')
ylim([1020.3, 1021.7])
lgd = legend(h(1:2), 'At 0.25 m', 'At 9 m');
set(lgd, 'location', 'southeast', 'fontsize', 12)


subplot(4, 2, 4) %Zonal geostrophic velocity
plot(middists./1000, ug(z(1),:),  'color', shallowCol, 'linestyle', '--', 'linewidth', 1)
plot(middists./1000, ug(z(2),:),  'color', deepCol, 'linestyle', '--', 'linewidth', 1)
ylabel('u_{geostrophic} [m/sec]')
% ylabel('geostrophic current [m/sec]')
ylim([-.1, .05])

subplot(4, 2, 5) %Horizontal buoyancy gradient
% plot(middists./1000, M2_x(1,:))
% hold on
% plot(middists./1000, M2_x(2,:))
hold on
plot(middists./1000, M2(z(1), :), 'color', shallowCol, 'linestyle', '--', 'linewidth', 1)
plot(middists./1000, M2(z(2), :), 'color', deepCol, 'linestyle', '--', 'linewidth', 1)
ylabel('M^2 [sec^{-2}]')

subplot(4, 2, 6) %Rossby number
plot(middists./1000, Ro(z(1), :),'color', shallowCol, 'linestyle', '--', 'linewidth', 1)%, 'Marker', 'o', 'MarkerFaceColor', shallowCol);
hold on
plot(middists./1000, Ro(z(2), :),  'color', deepCol, 'linestyle', '--', 'linewidth', 1)%, 'Marker', 'o', 'MarkerFaceColor', deepCol);
ylabel('Ro [\zeta / f]', 'fontsize', 14)

subplot(4, 2, 7) %Buoyancy frequency
plot(middists./1000, N2,'color', shallowCol, 'linestyle', '--', 'linewidth', 1)
hold on
ylabel('N^2 [sec^{-2}]')

subplot(4, 2, 8) %Richardson number
plot(middists./1000, Ri(1, :), 'color', shallowCol, 'linestyle', '--', 'linewidth', 1)
hold on
plot(middists./1000, Ri(2, :), 'color', deepCol, 'linestyle', '--', 'linewidth', 1)
ylabel('Ri [f^2 N^2 M^{-4}]')
ylim([0, 60])

for splot = 3:8
    subplot(4, 2, splot)
    xlabel('Distance [km]')
    xlim([0, 60])
    grid on; box on
end

% Plot section showing geeostrophic velocity calculated from Seaglider section
figure
[MIDDISTS,  ZMID] = meshgrid(middists, profiles.z(1:30));
[DISTS,  Z] = meshgrid(dists, profiles.z(1:30));
pcolor(MIDDISTS./1000, ZMID,ug)
shading interp
hold on
[C, H] = contour(DISTS./1000, Z, rho, [1020:.2:1028], 'color', 0.6 .* [1 1 1], 'linestyle', ':'); 
[C, H] = contour(DISTS./1000, Z, rho_smooth, [1020:.2:1028], 'color', 0.6 .* [1 1 1]); 
clabel(C, H, 'color', 0.6 .* [1 1 1])
ylim([0, intDepthInd])
set(gca,  'ydir', 'reverse');
cmocean('balance'); 
caxis([-.05, .05])
cb = colorbar;
ylabel(cb, 'u geostrophic [m/sec]', 'fontsize', 14)
scatter(dists./1000, 0.*ones(size(dists)), 'k', 'filled', 'v');  

