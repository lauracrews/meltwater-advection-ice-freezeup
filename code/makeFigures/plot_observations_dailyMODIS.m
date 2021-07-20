% Make plots comparing the in situ sea surface temperature observations 
% and the matched MODIS SST (matching was done by running
% match_observations_modis.m). Options to plot data from all 
%platforms on one figure or to make separate panels for Healy underway data, 
%Seaglider and uCTD data, and Wave Glider data. 

close all
clearvars -except AMSR2 profiles wvdata metData
plotPlatformComparison = false;

load('modisComparison.mat') %This was created above, slow
limits = [-1.7, 1]; %Axis limits

%Use this to select a time subset of the data for plotting
timePeriod = 0; %If timePeriod ~= 1, 2, 3, then will plot data from the entire time on one plot
if timePeriod == 0
    startTime = datenum('sept 19 2018')+.3; 
    endTime = datenum('oct 15 2018')-.4; 
elseif timePeriod == 1
    startTime = datenum('sept 19 2018'); 
    endTime = datenum('sept 24 2018');
    limits = [-1.7, 1]; %Different axis limits reflecting the differing min/max temperature values within each time period
elseif timePeriod == 2
    startTime = datenum('sept 30 2018'); 
    endTime = datenum('oct 7 2018'); 
    limits = [-1.5, 1];
elseif timePeriod == 3
    startTime = datenum('oct 8 2018'); 
    endTime = datenum('oct 15 2018'); 
    limits = [-2, 1.5];
end

inTimeMask = zeros(size(times));
inTimeMask(times >= startTime & times < endTime) = 1;

close all 
pedgs = -1.75:0.1:1.25; %Bin edges for creating the 2-d histograms
[X, Y] = meshgrid(pedgs, pedgs);

%Number of observations from each platform - use this to decide how to
%subsample
nHealy = length(find(platform == 1));
nProfile = length(find(platform == 2));
nWaveglider = length(find(platform == 3));

%Subset the Healy data
healyInc = round(nHealy / nProfile);
healySubsample = zeros(size(lons));
healySubsample(1:healyInc:nHealy) = 1; healySubsample(nHealy:end) = 1;

%Subset the Wave Glider data
wavegliderInc = 5;
wavegliderSubsample = zeros(size(lons));
wavegliderSubsample(1:(length(lons) - nWaveglider)) = 1; wavegliderSubsample((length(lons) - nWaveglider):wavegliderInc:end) = 1;

%Make a mask identifying which observations to use in the comparison
toPlotMask = zeros(size(lons));
toPlotMask(inTimeMask == 1 & healySubsample == 1 & wavegliderSubsample == 1) = 1;

%Make the 2d histogram comparing the in situ observations and MODIS SST
binCounts = hist3([temps(toPlotMask == 1), modisSST_nearObservations(toPlotMask == 1)], 'edges', {pedgs pedgs});
binCounts = binCounts'; %Need to transpose - see "Plot Histogram with Intensity Map" example in hist3 documentation
binCounts(binCounts == 0) = nan; %No color where the bins are empty

%Begin plotting
figure(4); set(gcf, 'color', 'w', 'pos', [176 534 944 414]) %Figure 4 - same as in plot_observation_SMOSsalinity.m to make side-by-side plots
subplot(1, 2, 1)
g = pcolor(X, Y, binCounts); set(g, 'linestyle', 'none')
hold on
colormap(gca, brewermap(20, 'greys'))

%Linear regression 
mdl = fitlm(temps(toPlotMask == 1), modisSST_nearObservations(toPlotMask == 1));
h = plot(mdl); 
delete(h(1)); %Remove the raw data points from the plot
set(h(2), 'color' ,'k', 'linewidth', 1) %Format linear fit line
set(h(3), 'color', 'k', 'linewidth', 1, 'linestyle', '--') %Format upper and lower confidence lines
set(h(4), 'color', 'k', 'linewidth', 1, 'linestyle', '--')
% title('Data from all platforms', 'fontsize', 12)
xlim(limits); ylim(limits);
x = limits(1):.1:limits(2);
y = x; plot(x, y, 'color', 0.7 .* [1 1 1], 'LineStyle', '-') %Plot 1:1 line
grid on
cb = colorbar;
ylabel(cb, 'Number of observations', 'fontsize', 14)
pbaspect([1 1 1])
ylabel(['MODIS sea surface temperature (', sprintf(char(176)), 'C)'], 'fontsize', 14)
xlabel(['Observed sea surface temperature (', sprintf(char(176)), 'C)'], 'fontsize', 14)
delete(legend) %This came with the plot(mdl) command
title('') %This also came with the plot(mdl) command
set(gca, 'fontsize', 12)

%% Make platform specific comparisons, if there are some observations in the time period
if plotPlatformComparison

    figure(1); set(gcf, 'color', 'w', 'pos', [188 563 1402 385])
    figure(2); set(gcf, 'color', 'w', 'pos', [192 166 1398 789])
    figure(3); set(gcf, 'color', 'w', 'pos', [188 563 1402 385])

    for curPlatform = 1:3

        if sum(platform == curPlatform & inTimeMask == 1) > 0 %There is some data to plot

            switch curPlatform
                case 1
                    cmap = 'Blues';
                    titleText = 'Healy underway data';
                case 2
                    cmap = 'Reds';
                    titleText = 'Seaglider and uCTD data';
                case 3
                    cmap = 'Greens';
                    titleText = 'Wave Glider data';
            end

            %Plot 2D histograms with linear fits
            figure(1); subplot(1, 3, curPlatform); hold on
            binCounts = hist3([temps(platform == curPlatform & inTimeMask == 1), modisSST_nearObservations(platform == curPlatform & inTimeMask == 1)], 'edges', {pedgs pedgs});
            binCounts = binCounts'; %Need to transpose - see "Plot Histogram with Intensity Map" example in hist3 documentation
            binCounts(binCounts == 0) = nan;
            g = pcolor(X, Y, binCounts); set(g, 'linestyle', 'none')
            colormap(gca, brewermap(20, cmap))
            mdl = fitlm(temps(platform == curPlatform & inTimeMask == 1), modisSST_nearObservations(platform == curPlatform & inTimeMask == 1));
            h = plot(mdl); delete(h(1));
            set(h(2), 'color' ,'k', 'linewidth', 1)
            set(h(3), 'color', 'k', 'linewidth', 1, 'linestyle', '--')
            set(h(4), 'color', 'k', 'linewidth', 1, 'linestyle', '--')
            title(titleText, 'fontsize', 12)

            figure(2); 
            %Bin matched MODIS data by day, plot histograms
            subplot(2, 3, curPlatform)
            edgs = min(modisTime_nearObservations(platform == curPlatform & inTimeMask == 1)):max(modisTime_nearObservations(platform == curPlatform & inTimeMask == 1));
            histogram(times(platform == curPlatform & inTimeMask == 1 & ~isnan(modisSST_nearObservations)), edgs);
            xlim([edgs(1) - 0.5, edgs(end)+.5])
            set(gca, 'XTick', edgs, 'XTickLabelRotation', 90)
            ylabel('Number of observations', 'fontsize', 12)
            title(titleText, 'fontsize', 12)

            %Bin in situ observations by day, plot histograms
            subplot(2, 3, curPlatform + 3)
            histogram(modisTime_nearObservations(platform == curPlatform & inTimeMask == 1 & ~isnan(modisSST_nearObservations)), edgs);
            xlim([edgs(1) - 0.5, edgs(end)+.5])
            set(gca, 'XTick', edgs, 'XTickLabelRotation', 90)
            ylabel('Number of matchhed MODIS SSTs', 'fontsize', 12)        

            %Bin observations by latitude, plot histograms
            figure(3); subplot(1, 3, curPlatform)
            edgs = round(min(lats(platform == curPlatform & inTimeMask == 1)), 1):0.2:round(max(lats(platform == curPlatform & inTimeMask == 1)), 1);
            histogram(lats(platform == curPlatform & inTimeMask == 1 & ~isnan(modisSST_nearObservations)), edgs);
            xlim([edgs(1) - 0.5, edgs(end)+.5])
            set(gca, 'XTick', edgs, 'XTickLabelRotation', 90)
            ylabel('Number of observations', 'fontsize', 12)
            title(titleText, 'fontsize', 12)
            xlabel('Latitude of observation', 'fontsize', 12)
        end
    end   

    %Format plots
    figure(1)
    for splot = 1:3
        subplot(1, 3, splot)
        hold on
        xlim(limits); ylim(limits);
        x = limits(1):.1:limits(2);
        y = x; plot(x, y, 'color', 0.7 .* [1 1 1], 'LineStyle', '-')
        grid on
        cb = colorbar;
        ylabel(cb, 'Number of observations', 'fontsize', 12)
        pbaspect([1 1 1])
        ylabel('MODIS SST', 'fontsize', 12)
        xlabel('Observed SST', 'fontsize', 12)
        delete(legend)
    end

    figure(2)
    for splot = 1:6
        subplot(2, 3, splot)
        datetick('x', 'mmm dd', 'keepticks', 'keeplimits')%, 'rotation', 90)
    end

end