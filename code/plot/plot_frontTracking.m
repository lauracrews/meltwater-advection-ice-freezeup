%Creates and plots simulated trajectories. Make vector timesereies of Ekman 
%velocities. Optionally animates the trajectories. 

% run plot_DOT_geostrophic first to calculate geostrophic velocities on a 
% regularly spaced lat/lon grid and plot DOT and geostrophic velocity vectors

close all; 
clearvars -except rootPath AMSR2 profiles wvdata metData ug vg lon_sub lat_sub

saveFigs = false;
saveDir = [rootPath, 'figures/fig7/'];

% Convert geostrophic velocities (calculated in plot_DOT_geostrophic.m) to column vectors 
ug = ug(:); vg = vg(:);
%%
makeAnimation = false; 
[~, colors] = defineSODAconstants;
rhoWater = 1026; %Reference density of Seawatere kg/m^3
D_ek = 10; %Assumed Ekman depth

%Make sure projection matches that of the figure on which the tracks will be plotted
%(i.e. the map made by plot_DOT_geostrophic.m) 
minlon = -152; maxlon = -142; minlat = 72; maxlat = 75.7;
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

%Starting locations for trajectories
startPos = [-147.125, 73; -147.125, 73.05; -147.125, 73.1; -147.25, 73.15; -147.375, 73.15; -147.5, 73.15; -147.625, 73.15; -147.75, 73.15; -147.875, 73.2; -148, 73.2];

%Time period for trajectories
startTime = datenum('sept 21, 2018');
endTime = datenum('oct 15 2018')-(1/24); %stop on the last timestep of the previous day

%Load atmospheric forcing, for calculating Ekman velocties
filename = 'ERA5.nc';
eraTime = double(ncread(filename, 'time'))/24 + datenum('1900-01-01 00:00:00');    
[~, startInd] = min(abs(eraTime - startTime));
[~, endInd] = min(abs(eraTime - endTime)); 
trajTimesteps = eraTime(startInd:endInd);
dt = mean(diff(trajTimesteps)) * 60 * 60 * 24; %convert timesteps (1/24 days) to seconds

eraLon = double(ncread(filename, 'longitude'));
eraLat = double(ncread(filename, 'latitude'));

[~, startLon] = min(abs(eraLon - minlon));
[~, endLon] = min(abs(eraLon - maxlon));
[~, startLat] = min(abs(eraLat - maxlat));
[~, endLat] = min(abs(eraLat - minlat));

eraLon = double(ncread(filename, 'longitude', startLon, [endLon - startLon + 1]));
eraLat = double(ncread(filename, 'latitude', startLat, [endLat - startLat + 1]));
[eraLon, eraLat] = meshgrid(eraLon, eraLat);
eraLon = eraLon'; eraLat = eraLat';

%Load ERA5 meridional and zonal turbulent surface stress
metss = double(ncread(filename, 'metss', [startLon startLat startInd],...
    [endLon - startLon + 1,  endLat - startLat + 1, endInd - startInd + 1]));

mntss = double(ncread(filename, 'mntss', [startLon startLat startInd],...
    [endLon - startLon + 1,  endLat - startLat + 1, endInd - startInd + 1]));

%% Filter the surface stress data to retain frequencies of at least 1 day
dt_f = mean(diff(eraTime)) * 24; %Sampling interval (hours). Need to distinguish from dt, which is the timesep in seconds for the model trajectories
Fs = 1/dt_f; %Sampling frequency
fc = 1/16; %Cutoff frequency, 1/Hours
Wn = fc/(Fs/2);
ord = 5; %Define filter order
[b,a] = butter(ord, Wn);

mntss_f = nan(size(mntss));
metss_f = mntss_f;
[I, J, K] = size(metss);
for i = 1:I  
    for j = 1:J
        metss_f(i, j, :) = filtfilt(b, a, squeeze(metss(i, j, :)));
        mntss_f(i, j, :) = filtfilt(b, a, squeeze(mntss(i, j, :)));  
    end
end

%%
%Empty arrays to hold the tracking results

%Trajectories for combined velocity
trajLons_tot = nan.*ones(length(trajTimesteps), size(startPos, 1));
trajLats_tot = nan.*ones(length(trajTimesteps), size(startPos, 1));

%Trajectories for geostrophic velocity only
trajLons_g = nan.*ones(length(trajTimesteps), size(startPos, 1));
trajLats_g = nan.*ones(length(trajTimesteps), size(startPos, 1));

%Trajectories for Ekman velocity only
trajLons_e = nan.*ones(length(trajTimesteps), size(startPos, 1));
trajLats_e = nan.*ones(length(trajTimesteps), size(startPos, 1));

%Ekman velocities
allue = nan.*ones(length(trajTimesteps), size(startPos, 1)); 
allve = nan.*ones(length(trajTimesteps), size(startPos, 1));

%Surface stresses
alltaux = nan.*ones(length(trajTimesteps), size(startPos, 1));
alltauy = nan.*ones(length(trajTimesteps), size(startPos, 1));
%%
for startNum = 1:size(startPos, 1) %Create multiple trajectories - one for each start position
    disp(['Calculating simulated trajectory ', num2str(startNum), ' of ', num2str(size(startPos, 1))])
    
    for vel = 1 %Advect by 1) combined velocity, 2) geostrophic velocity only, or 3) Ekman velocity only

        clear u_cur v_cur ug_cur vg_cur ue_cur v_cur
        curLon = startPos(startNum, 1); curLat = startPos(startNum, 2); %Initialize tracking at this location

        for curTime = 1:length(trajTimesteps) %Iterate through all the timesteps

            if vel ~= 3 %No need to calculate geostrophic velocity if advecting with only Ekman velocity

                %Find the closest geostrophic velocity point. lon_sub and lat_sub were
                %made in plot_DOT_geostrophic and are the coordinates of the
                %geostrophic velocity data points
                dists  = zeros(size(lon_sub));
                for i = 1:size(lon_sub, 1)
                    for j = 1:size(lon_sub, 2)
                        dists(i, j) = m_lldist([curLon, lon_sub(i, j)], [curLat, lat_sub(i, j)]);
                    end
                end
                [~, ind] = min(dists(:));
                ug_cur = ug(ind); vg_cur = vg(ind); %Geostrophic velocity at the current time
                clear dists
            end

            if vel ~= 2 %No need to calculate Ekman velocity if advecting with only geostrophic velocity

                %ERA5 surface stress data at the current time step in the whole region 
                metss_f_cur = squeeze(metss_f(:, :, curTime)); 
                mntss_f_cur = squeeze(mntss_f(:, :, curTime)); 

                %Identify the spatially-closest ERA5 data
                dists  = nan .* ones(size(eraLon));
                for i = 1:size(eraLon, 1)
                    for j = 1:size(eraLon, 2)
                        dists(i, j) = m_lldist([curLon, eraLon(i, j)], [curLat, eraLat(i, j)]);
                    end
                end
                [~, ind] = min(abs(dists(:))); clear dists
                tau_airwater_x = metss_f_cur(ind); tau_airwater_y = mntss_f_cur(ind);

                %Calculate average ocean velocity in the Ekman layer
                ve_cur =  -(1/gsw_f(curLat)) * tau_airwater_x * (1/(rhoWater * D_ek));
                ue_cur = (1/gsw_f(curLat)) * tau_airwater_y * (1/(rhoWater * D_ek));
            end

            switch vel %Determine the velocity at this timestep
                case 1 %Combined velocity
                    u_cur = ug_cur + ue_cur;
                    v_cur = vg_cur + ve_cur;       
                case 2 %Gestrophic only
                    u_cur = ug_cur;
                    v_cur = vg_cur;
                case 3 %Ekman only
                    u_cur = ue_cur;
                    v_cur = ve_cur;
            end

            %Calculate new position
            displacement = dt * hypot(u_cur, v_cur); %Speed [m/sec] * time [sec]
            heading = calculateCompassAngle(u_cur, v_cur); %Compass heading of the velocity, needed for m_fdist
            [curLon, curLat] = m_fdist(curLon, curLat, heading, displacement); curLon = curLon - 360;

            switch vel %Update arrays
                case 1
                    trajLons_tot(curTime, startNum) = curLon;
                    trajLats_tot(curTime, startNum) = curLat;
                case 2
                    trajLons_g(curTime, startNum) = curLon;
                    trajLats_g(curTime, startNum) = curLat;
                case 3
                    trajLons_e(curTime, startNum) = curLon;
                    trajLats_e(curTime, startNum) = curLat;

                    allue(curTime, startNum) = ue_cur;
                    allve(curTime, startNum) = ve_cur;

                    alltaux(curTime, startNum) = tau_airwater_x;
                    alltauy(curTime, startNum) = tau_airwater_y;
            end

            %Optional plotting code to make sure code is calculating
            %everything correctly 
    %         plotEvery = 1;
    %         eraLon = eraLon(:); eraLat = eraLat(:);
    %         q(1) = m_quiver(eraLon(1:plotEvery:end), eraLat(1:plotEvery:end), metss_cur(1:plotEvery:end), mntss_cur(1:plotEvery:end), 0, 'k');
    %         hold on
    %         q(2) = m_quiver(curLon, curLat, alltaux(curTime, :), alltauy(curTime, :), 0, 'g', 'linewidth', 1);
    %         q(3) = m_quiver(curLon, curLat, allue(curTime, :), allve(curTime, :), 0, 'r', 'linewidth', 1);
    % 
    %         h(3) = m_scatter(driftLons_e(:), driftLats_e(:), 12, 'k', 'filled');
    %         title(datestr(driftTimes(curTime), 'dd mmm hh:MM'))
    %         m_grid
    %         drawnow
    %         clf 
        end
    end
end

%% Optionally open a different figure to plot the traajectories on
% openfig('ModisTerraNSST_sept30toOct7_withIce.fig')
% openfig('ModisTerraNSST_sept14to21_withIce_wavegliders.fig')
openfig([rootPath, 'figures/fig7/DOT_geostrophic.fig'])
% openfig('frontTracking_ekmanDepth15.fig')

%Plot trajectories, with colors changing at important times
days = [datenum('sept 26 2018'), datenum('oct 3 2018'), datenum('oct 11 2018'), endTime]; 
prevInd = 1;
clrs = [0 0 0; 0.3 0.3 0.3; 0.6 0.6 0.6; 0.9 0.9 0.9]; %Array of colors (grayscale) to plot the trajectories
for day = 1:length(days)
    [~, ind] = min(abs(trajTimesteps - days(day)));
    lons = trajLons_tot(prevInd:ind, :); lats = trajLats_tot(prevInd:ind, :); %Choose which trajectory (Ekman, geostrophic, combined) by editing here
    h(3) = m_scatter(lons(:), lats(:), 10, clrs(day, :), 'filled');
%     m_text(lons(end, 1), lats(end, 1), ['   ', datestr(driftTimes(ind), 'dd mmm')], 'fontsize', 12)
    prevInd = ind;
end

%Add trajectory start positions to the map
m_scatter(startPos(:, 1), startPos(:, 2), 20, 'w', 'filled');
m_text(startPos(1, 1), startPos(1, 2), ['  ', datestr(startTime(1), 'dd mmm')], 'fontsize', 12, 'color', 'w')

 
%% Plot crossings of the 26 g/kg isohaline - code copied from plot_surfaceSalinity_byTime.m
saltOutcrop = 26;

%Load uCTD and Seaglider profiles
% profiles = loadProfiles;

%Load underway data from the Healy system 
% metData = load_correct_underwayData;

%Identify 26 salinity outcrop = front location in ship data
startTime_shipFront = datenum('sept 26 2018') + .4;
endTime_shipFront = datenum('sept 27 2018') - .4;
startInd = find(abs(metData.times - startTime_shipFront) == min(abs(metData.times - startTime_shipFront)));
endInd = find(abs(metData.times - endTime_shipFront) == min(abs(metData.times - endTime_shipFront)));

times_met =  metData.times(startInd:endInd); 
lats_met = metData.lats(startInd:endInd); lons_met = metData.lons(startInd:endInd); 

ind = find(abs(metData.SAs(startInd:endInd) - saltOutcrop) == min(abs(metData.SAs(startInd:endInd) - saltOutcrop)));
outcropPos(1, :) = [lons_met(ind), lats_met(ind)]; 
outcropTime(1) = times_met(ind);

%Add front location from ship to map
% m_scatter(outcropPos(1, 1), outcropPos(1, 2), 150, '^', 'w', 'filled')
% m_scatter(outcropPos(1, 1), outcropPos(1, 2), 100, clrs(2,:), '^', 'filled')
% m_text(outcropPos(1, 1) + .2, outcropPos(1, 2), ['  ', datestr(outcropTime(1), 'dd mmm')], 'color', 'w', 'fontsize', 12)
  %%        
for transect = 2:3
    switch transect
        case 2
            %Identify front location on the Seaglider 199 crossing of the front
            startTime = datenum('sept 30 2018'); endTime = datenum('oct 6 2018');  
        case 3
            startTime = datenum('oct 6 2018'); endTime = datenum('oct 11 2018') + .2; %Should match endTime in plot_surfaceSalinity_byTime transect 3
    end
    
    profNums = find(strcmp(profiles.dataset, 'sg199') & profiles.times >= startTime & profiles.times <= endTime);
    sg_lons = profiles.lons(profNums); sg_lats = profiles.lats(profNums);
    times = profiles.times(profNums);
    salts = nanmean(profiles.SA(1:5, profNums));
        
     %Identify 26 salinity outcrop = front location in Seaglider data
    [~, indFront] = min(abs(salts - saltOutcrop));

    if isempty(indFront) %This should not occur, though the isohaline barely outcrops on Transect 3
        indFront = length(salts); %Defaults to the end of the transect
    end

    outcropPos(transect, :) = [sg_lons(indFront), sg_lats(indFront)]; 
    outcropTime(transect) = times(indFront);

end

m_line([startPos(4, 1), outcropPos(1, 1)], [startPos(4, 2), outcropPos(1, 2)], 'linewidth', 4, 'linestyle', ':', 'color', clrs(1, :))
m_line(outcropPos(1:2, 1), outcropPos(1:2, 2), 'linewidth', 4, 'linestyle', ':', 'color', clrs(2, :))
m_line(outcropPos(2:3, 1), outcropPos(2:3, 2), 'linewidth', 4, 'linestyle', ':', 'color', clrs(3, :))

for i = 1:3
    m_scatter(outcropPos(i, 1), outcropPos(i, 2), 150, 'w', '^', 'filled')
    m_scatter(outcropPos(i, 1), outcropPos(i, 2), 100, clrs(i+1,:), '^', 'filled')
    m_text(outcropPos(i, 1) + .2, outcropPos(i, 2), ['  ', datestr(outcropTime(i), 'dd mmm')], 'color', 'w', 'fontsize', 12)
end

transectEndTime = times(end);

%%
%Save map figure with simulated trajectories
saveName = 'frontTracking';
if saveFigs
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end      

%% Average wind and Ekman velocity in a box, for plotting vector timeseries

%Region in which to average velocities
pt1 = [-149, 73]; pt2 = [-145, 73]; pt3 = [-145, 74.5]; pt4 = [-149, 74.5]; 
% pt1 = [-147, 74.5]; pt2 = [-145, 74.5]; pt3 = [-145, 75.5]; pt4 = [-147, 75.5]; 

pts = [pt1(1), pt2(1), pt3(1), pt4(1); pt1(2), pt2(2), pt3(2), pt4(2)];

%Plot box of area used for averaging wind - need to select map figure
% figure(1); subplot(6, 1, 1:4, ax2)
[in, on] = inpolygon(eraLon, eraLat, pts(1, :), pts(2, :));
% m_contour(eraLon, eraLat, double(in), [0, 1], 'color', 'k', 'linewidth', 2, 'linestyle', '--')

startTime = datenum('sept 21, 2018');
endTime = datenum('oct 13 2018')-.125; %('oct 15 2018')-.125; %stop on the last timestep of the previous day

%Wind and surface stress from ERA5
atmFluxes = makeERA5timeseries(pts, startTime, endTime, {'u10', 'v10', 'metss', 'mntss'});
tau_airwater_x = [atmFluxes.metss]'; tau_airwater_y = [atmFluxes.mntss]';
times = [(atmFluxes.time)]';
windSpeeds = sqrt([atmFluxes.u10].^2 + [atmFluxes.u10].^2);

dt = diff([atmFluxes(1:2).time]) * 24; %Sampling interval (hours)
Fs = 1/dt; %Sampling frequency
fc = 1/24; %Cutoff frequency, 1/Hours
Wn = fc/(Fs/2);
ord = 5;
[b,a] = butter(ord, Wn);

u10low = filtfilt(b, a, [atmFluxes.u10]);
v10low = filtfilt(b, a, [atmFluxes.v10]);
tau_airwater_x = filtfilt(b, a, tau_airwater_x);
tau_airwater_y = filtfilt(b, a, tau_airwater_y);

%Calculate Ekman velocitties
v_ek = -(1/gsw_f(mean(pts(2, :)))) * tau_airwater_x * (1/(rhoWater * D_ek));
u_ek = (1/gsw_f(mean(pts(2, :)))) * tau_airwater_y * (1/(rhoWater * D_ek));

% Plot wind and Ekman drift vectors
% figure; set(gcf, 'color', 'w', 'pos', [48 92 1492 863])
subplot(6, 1, 5:6); hold on
% subplot(3, 1, 1); hold on
pt1 = quiver(times, zeros(size(times)), u10low' ./ 10, v10low' ./ 10, 0, 'linewidth', 1, 'color', colors.blue, 'ShowArrowHead', 'off');
set(pt1, 'autoscale', 'off'); axis equal 
pt2 = quiver(datenum('Oct 1 2018'), -0.7, 10 ./ 10, 0 ./ 10, 0, 'linewidth', 1, 'color', colors.blue);
set(pt2, 'autoscale', 'off'); axis equal
text(datenum('Oct 1 2018'), -0.5, '10-m wind, 10 m/sec', 'fontsize', 12, 'color', colors.blue);
pt3 = quiver(times, zeros(size(times)), u_ek .* 10, v_ek .* 10, 0, 'linewidth', 1, 'color', colors.orange, 'ShowArrowHead', 'off'); 
set(pt3, 'autoscale', 'off'); axis equal 
pt4 = quiver(datenum('Oct 4 2018'), -0.7, 1, 0, 0, 'linewidth', 1, 'color', colors.orange);
set(pt2, 'autoscale', 'off'); axis equal 
text(datenum('Oct 4 2018'), -0.5, 'Mean Ekman velocity, 10 cm/sec', 'fontsize', 12, 'color', colors.orange);

xlim([min(times) - .5, max(times) + .5])
ylim([-2, 2])
grid on; box on
set(gca, 'ytick', '', 'xtick', startTime:1:ceil(endTime), 'fontsize', 12, 'XTickLabelRotation', 90)
datetick('x', 'mmm dd', 'keeplimits', 'keepticks')

%Make a compass rose in the upper left corner of plot
line(datenum('sep 21 2018') + [.1, 0.9], 1.3 .* [1 1], 'color', 'k', 'linewidth', 1)
line(datenum('sep 21 2018') + [.5, .5], [.95 1.65], 'color', 'k', 'linewidth', 1)
text(datenum('sep 21 2018') + .4, 1.82, 'N')
text(datenum('sep 21 2018') + .4, .78, 'S')
text(datenum('sep 21 2018') - .22, 1.3, 'W')
text(datenum('sep 21 2018') + .93, 1.3, 'E')
%%
%Additional plots showing one component direction wind and Ekman velocities 
figure; set(gcf, 'color', 'w', 'pos', [48 92 1492 863])

subplot(3, 1, 2)
yyaxis left; hold on
plot(times, [atmFluxes.u10]')
plot(times, u10low')

yyaxis right
plot(times, -v_ek)%, 'color', 'r')

subplot(3, 1, 3)
yyaxis left; hold on
plot(times, [atmFluxes.v10]')
plot(times, v10low')

yyaxis right; hold on
plot(times, u_ek)%, 'color', 'r')

for splot = 2:3 %Format plots
    subplot(3, 1, splot)

    yyaxis left
    ylim([-12, 12])
    ylabel('Wind velocity [m/sec]')
    
    yyaxis right; hold on
    ylim(.15 .* [-1, 1])
    ylabel('Ekman velocity [m/sec]')
    
    xlim([min(times) - .5, max(times) + .5])
    datetick('x', 'mm/dd', 'keeplimits')
    legend('v wind at 10 m', 'u Ekman')
    grid on
end

%Adjust vector plot position to align with lower plots
ps2 = get(gca, 'pos');

subplot(3, 1, 1)
ps1 = get(gca, 'pos');
ps1(1) = ps2(1);
ps1(3) = ps2(3);
set(gca, 'pos', ps1)

%Save figure with wind and Ekman transport vectors
saveName = 'ekmanVectors';
if saveFigs    
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end

%% Make an animation of drift trajectories
if makeAnimation
    v1 = VideoWriter([rootPath, 'figures/fig7/frontTracking_combined_ekmanDepth10m_test.mp4'], 'MPEG-4'); %
    v1.Quality = 100; v1.FrameRate = 36; %frames / sec 
    open(v1)
    
    trajLons_toPlot = trajLons_tot; %Edit acording to which trajectories you want to plot
    trajLats_toPlot = trajLats_tot;
    
    plot_surfaceSalinity_byTime; %Run this to identify the meltwater front in observations
    %%% This ^ alters the map projection, be sure to update! 
    minlon = -152; maxlon = -142; minlat = 72; maxlat = 75.7; %Make sure projection matches that of the figure on which the tracks will be plotted
    m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

    %Load sea ice data to add to animation
%     load AMSR2_2018.mat
    amsr_lon = AMSR2.lon; amsr_lat = AMSR2.lat; 

    %Start with the earlier sea surface temperature image
    close all
    openfig('ModisTerraNSST_sept14to21_wavegliders.fig')

    %Optional code to make sure the lat/lon subsetting when loading is
    %operating correctly. If it is, these stress vectors should exactly match
    %the metss and mentss vectors
    % eraLon_test = double(ncread(filename, 'longitude'));
    % eraLat_test = double(ncread(filename, 'latitude'));
    % 
    % [eraLon_test, eraLat_test] = meshgrid(eraLon_test, eraLat_test);
    % eraLon_test = eraLon_test'; eraLat_test = eraLat_test';
    % 
    % metss_test = double(ncread(filename, 'metss', [1 1 startInd],...
    %     [Inf,  Inf, endInd - startInd + 1]));
    % 
    % mntss_test = double(ncread(filename, 'mntss', [1 1 startInd],...
    %     [Inf,  Inf, endInd - startInd + 1]));

    
    %Create color vectors that change at the important times
    days = [datenum('sept 26 2018'), datenum('oct 3 2018'), datenum('oct 11 2018'), endTime]; 
    prevInd = 1;
    clrs = [0 0 0; 0.3 0.3 0.3; 0.6 0.6 0.6; 0.9 0.9 0.9]; %0.5 0.5 0.5;
    colVec = zeros(length(trajTimesteps), 3);
    for day = 1:length(days)
        [~, ind] = min(abs(trajTimesteps - days(day)));
        colVec(prevInd:ind, :) = repmat(clrs(day, :), length(prevInd:ind), 1);
        prevInd = ind;
    end

    %Iterate through timesteps and plot that frame
    for curTime = 1:length(trajTimesteps) 

        %Change to the later MODIS SST image
        if trajTimesteps(curTime) == datenum('Oct 7 2018')
            close all
            openfig('ModisTerraNSST_sept30toOct7.fig')
        end

        %Add today's sea ice concentration to the plot
        [~, iceInd] = min(abs(AMSR2.mattime - floor(trajTimesteps(curTime))));
        curIce = AMSR2.SIC(:, :, iceInd);
        curIce(curIce < 0.15) = nan;
        h(4) = m_pcolor(AMSR2.lon, AMSR2.lat, curIce);
        set(h(4), 'facealpha', .9)
        cmocean('ice'); caxis([0, 1])
        hold on
        title(datestr(trajTimesteps(curTime), 'dd mmm hh:MM'))
        
        %Add trajectories up to this time 
        curLons = trajLons_toPlot(1:curTime, :); curLats = trajLats_toPlot(1:curTime, :);
        h(3) = m_scatter(curLons(:), curLats(:), 12, repmat(colVec(1:curTime, :), 10, 1), 'filled');
        h(2) = m_scatter(curLons(end, :), curLats(end, :), 15, 'w', 'filled');

        %Optionally add vectors of surface stress to make sure code is
        %working correctly
%         metss_cur = metss(:, :, curTime); mntss_cur = mntss(:, :, curTime); 
%         plotEvery = 1; %Interval at which to plot stress vectors
    %     q(1) = m_quiver(eraLon(1:plotEvery:end), eraLat(1:plotEvery:end), metss_cur(1:plotEvery:end), mntss_cur(1:plotEvery:end), 0, 'k', 'linewidth', 1);
    %     q(3) = m_quiver(driftLons_toPlot(curTime, :), driftLats_toPlot(curTime, :), alltaux(curTime, :), alltauy(curTime, :), 0, 'g', 'linewidth', 1);
    %     q(4) = m_quiver(driftLons_toPlot(curTime, :), driftLats_toPlot(curTime, :), allue(curTime, :), allve(curTime, :), 0, 'r', 'linewidth', 1);

        %Add observed locations of the meltwater edge, if they have
        %occurred by this timestep
        if trajTimesteps(curTime) >= outcropTime_ship
                %Add front location from ship to map
            h(5) = m_scatter(outcrop_ship(1), outcrop_ship(2), 150, '^', 'w', 'filled');
            h(6) = m_scatter(outcrop_ship(1), outcrop_ship(2), 100, clrs(2,:), '^', 'filled');
            h(7) = m_text(outcrop_ship(1) + 1, outcrop_ship(2), ['  Healy, ', datestr(outcropTime_ship, 'dd mmm')], 'color', 'k', 'fontsize', 12);
        end
        
        if trajTimesteps(curTime) >= outcropTime_sg
            h(9) = m_scatter(outcrop_sg(1), outcrop_sg(2), 150, 'w', '^', 'filled');
            h(10) = m_scatter(outcrop_sg(1), outcrop_sg(2), 100, clrs(3,:), '^', 'filled');
            h(11) = m_text(outcrop_sg(1) + 1, outcrop_sg(2), ['  SG199, ', datestr(outcropTime_sg, 'dd mmm')], 'color', 'k', 'fontsize', 12);
        end
        
        if trajTimesteps(curTime) >= transectEndTime
            h(12) = m_scatter(sg_lons(end), sg_lats(end), 150, 'w', '^', 'filled');
            h(13) = m_scatter(sg_lons(end), sg_lats(end), 100, clrs(4,:), '^', 'filled');
            h(14) = m_text(-146, 74.25, ['  End of SG199 transect', sprintf('\n'), '   through front, ', datestr(transectEndTime, 'dd mmm')], 'color', 'k', 'fontsize', 12);
        end

        F1 = getframe(gcf);
        writeVideo(v1, F1) 
        delete(h)
    %     delete(q)
    end
    close(v1);
end
