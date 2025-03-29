%Plot map of melt-out and freeze-up date from AMSR2
%Adapted from code from Luc Rainville

close all;

iceThresh = 0.15;
dayThreshold = 5; %Number of days in a fow that that a cell must be ice covered or ice free
[moorings, colors] = defineSODAconstants; %Mooring locations and standard colors
% load('Beaufort_AMSR2_2018.mat') %Load AMSR2 data
sz = size(AMSR2.lon);

saveFigs = false;
saveDir = [rootPath ,'figures/fig1/'];
saveName = 'ice_formation_melt';

%% Calculate days of melt out and freeze up

% first day when the ice concentration is >0.15 for more than DT days.
meltOut = nan(size(AMSR2.lon)); %length(bremen.lat), length(bremen.lon))*NaN;

% last day when the ice concentration is >0.15 for more than DT days.
freezeUp = nan(size(AMSR2.lon)); %zeros(length(bremen.lat), length(bremen.lon))*NaN;

alwaysOpen = zeros(size(AMSR2.lon)); %zeros(length(bremen.lat), length(bremen.lon));
alwaysClosed = zeros(size(AMSR2.lon));

for k1=1:sz(1)
    for k2=1:sz(2)
        iceTimeseries = AMSR2.SIC(k1, k2, :); iceTimeseries = iceTimeseries(:); %iceTimeseries(1:startInd) = nan;
        if sum(~isnan(iceTimeseries) == 0)
            continue
        end
        meltind = find(iceTimeseries < iceThresh);%, 1, 'first'); %Freeze-up day in each grid cell
        freezeind = find(iceTimeseries >= iceThresh);%, 1, 'first'); %Freeze-up day in each grid cell
        
        if all(iceTimeseries < iceThresh)
            alwaysOpen(k1, k2) = 1;
            continue
        elseif all(iceTimeseries >= iceThresh)
            alwaysClosed(k1, k2) = 1;
            continue
        end
        
        %Work backward to find the last day
        if length(meltind) >= dayThreshold
            for i = 1:length(meltind)-dayThreshold+1
                if all(diff(meltind(i:i+dayThreshold-1)) == 1)% & sum(diff(meltind(i:i+dayThreshold-1))) == dayThreshold %all(diff(meltind(i-dayThreshold+1:i)) == 1)
                    meltOut(k1, k2) = AMSR2.mattime(meltind(i));
                    freezeind(freezeind < meltind(i)) = [];
                    break
                end
            end
        end
        
        if isnan(meltOut(k1, k2))
            alwaysClosed(k1, k2) = 1;
            continue
        end
        
        freezeind(AMSR2.mattime(freezeind) < datenum('sept 21, 2018')) = [];
        if length(freezeind) >= dayThreshold
            for i = 1:length(freezeind)-dayThreshold+1
                if all(diff(freezeind(i:i+dayThreshold-1)) == 1) %& sum(diff(freezeind(i:i+dayThreshold-1))) == dayThreshold
                    freezeUp(k1, k2) = AMSR2.mattime(freezeind(i));
                    break
                end
            end
        end
        
        if isnan(freezeUp(k1, k2))
            alwaysOpen(k1, k2) = 1;
            continue
        end
        
        
        makePlot = false;
        if makePlot & alwaysOpen(k1, k2) == 0 && alwaysClosed(k1, k2)==0
            close all
            plot(AMSR2.mattime, iceTimeseries);
            hold on
            line(meltOut(k1, k2) .* [1 1], [0, 1], 'color', 'k', 'linewidth', 1);
            line(freezeUp(k1, k2) .* [1 1], [0, 1], 'color', 'r', 'linewidth', 1);
        end
     
    end
end
%%
%Calculate first and last day (or detect if always ice covered) for each
%lat/lon grid point
% for k1=1:sz(1)
%     for k2=1:sz(2)
%         ii = conv2P(double(squeeze(AMSR2.SIC(k1,k2,:))<0.15), ones(dayThreshold,1)/dayThreshold); % this is one only if it was true DT/2 days before and after
%         if ~isempty(find(ii==1, 1))
%             meltOut(k1,k2)=AMSR2.mattime(find(ii==1,1));
%             freezeUp(k1,k2)=AMSR2.mattime(find(ii==1,1,'last'));
%         else 
%             minIce(k1,k2)=1; % completely ice covered the whole time. 
%         end
%     end
% end
% 
% 
% 
% 
% %Iterate through all AMSR2 grid cells
% firstIce = nan .* ones(size(AMSR2.lat));
% for i = size(AMSR2.lat, 1):-1:1
%     for j = size(AMSR2.lat, 2):-1:1      
%         iceTimeseries = AMSR2.SIC(i, j, :); iceTimeseries = iceTimeseries(:); iceTimeseries(1:startInd) = nan;
%         ind = find(iceTimeseries >= iceThresh, 1, 'first'); %Freeze-up day in each grid cell
%         if ~isempty(ind); firstIce(i, j) = AMSR2.mattime(ind); end
%     end
% end


%%
close all

lonlimits = [-165 -135];
latlimits = [69.5 78];
m_proj('Lambert','lon',lonlimits,'lat',latlimits)
year = 2018;

%Limits for colormaps
colorIncrement_meltout = 7; %step size for colormap, in days
colorLims_meltOut = [datenum(year,7,15) datenum(year,9,15)]; %Earliest and latest day for color axis

colorIncrement_freezeUp = 3;%2; %step size for colormap, in days
colorLims_freezeUp = [datenum(year,10,9) datenum(year,10,22)];
colorLims_freezeUp = [datenum(year,10,11) datenum(year,10,27)];

%Begin plotting melt out timing
meltOut_plot = meltOut;
meltOut_plot(alwaysOpen == 1) = colorLims_meltOut(1);
meltOut_plot(meltOut < colorLims_meltOut(1)) = colorLims_meltOut(1);
meltOut_plot(alwaysClosed == 1) = nan;

figure(1); set(gcf, 'color', 'w'); paper(16,7)
h1=subplot_m(1, 2, 1, 0.05);
m_pcolor(AMSR2.lon, AMSR2.lat, meltOut_plot); shading flat;
caxis([colorLims_meltOut(1)-colorIncrement_meltout/2 colorLims_meltOut(2)+colorIncrement_meltout/2]);
title({'Date of ice melt', num2str(year)}, 'fontsize',20)

%Set up melt-out colormap
c_tick = colorLims_meltOut(1):colorIncrement_meltout:colorLims_meltOut(2)+colorIncrement_meltout/2;
colorLims_meltOut=[colorLims_meltOut(1) c_tick(end)];
cmap = flipud(cbrewer('seq','YlGnBu', length(c_tick)));

%Make colorbar
c1 = colorbar_out(0.5, 0, 0.02);
c_tick = colorLims_meltOut(1):colorIncrement_meltout:colorLims_meltOut(2)+colorIncrement_meltout/2;
clear c_tick_label; count = 1;
for k=1:2:length(c_tick)
  c_tick_label{count} = datestr(c_tick(k), 'mmm dd'); %Create vector of colorbar labels
  count = count + 1;
end
pp =get(c1,'position');
set(c1, 'ytick',c_tick(1:2:end), 'YTickLabel', c_tick_label, 'position',[pp(1)-pp(3),pp(2)+pp(4)/3,pp(3:4)]);
colormap(h1,cmap)

%Begin plotting freeze up timing
freezeUp_plot = freezeUp;
freezeUp_plot(alwaysOpen == 1) = colorLims_freezeUp(2);
freezeUp_plot(freezeUp<colorLims_freezeUp(1))=colorLims_freezeUp(1);
freezeUp_plot(alwaysClosed == 1) = nan;

h2=subplot_m(1, 2, 2, 0.05);
m_pcolor(AMSR2.lon, AMSR2.lat, freezeUp_plot); shading flat;
caxis([colorLims_freezeUp(1)-colorIncrement_freezeUp/2 colorLims_freezeUp(2)+colorIncrement_freezeUp/2]);
caxis([colorLims_freezeUp(1) - 0.5*colorIncrement_freezeUp, colorLims_freezeUp(2) + 0.5*colorIncrement_freezeUp]);

title({'Date of ice formation',num2str(year)},'fontsize',20)

%Set up freeze up colormap
c_tick = colorLims_freezeUp(1):colorIncrement_freezeUp:colorLims_freezeUp(2)+colorIncrement_freezeUp/2;
colorLims_freezeUp=[colorLims_freezeUp(1) c_tick(end)];
cmap2 = cbrewer('seq','Greys', length(c_tick)); %

%Make colorbar for freeze up time
c2 = colorbar_out(0.5,0, 0.02);
clear c_tick_label; count = 1;
for k=1:2:length(c_tick) %Create vector of colorbar labels, labelling every other color increment
  c_tick_label{count} = datestr(c_tick(k),'mmm dd');
  count = count + 1;
end
pp = get(c2, 'position');
set(c2,'ytick', c_tick(1:2:end), 'YTickLabel',c_tick_label, 'position',[pp(1)-pp(3),pp(2)+pp(4)/3,pp(3:4)]);
colormap(h2,cmap2)
%%
%Additions to both plots
for splot = 1:2
    h2=subplot_m(1, 2, splot, 0.05);
    hold on
    m_contour(ibcao.LON, ibcao.LAT,ibcao.d,[-1000 -1000],'color',[1 1 1]*0.9,'linewidth',1); %Plot bathymetry
    m_contour(AMSR2.lon, AMSR2.lat, alwaysClosed,[1 1]*0.15,'color',[0 0 1],'linewidth',1)
    m_gshhs_i('patch', 0.7*[1 1 1]); %Draw coastline

    m_grid('xtick', [-180:10:-120], 'tickdir', 'out', 'ytick', [60:5:85], 'linest', '-', 'fontsize', 12);
    
    %Add box for temperature and salinity plots
    minlon_box = -150; maxlon_box = -143; minlat_box = 72; maxlat_box = 75.5;
    pt1 = [minlon_box, minlat_box]; pt2 = [maxlon_box, minlat_box]; pt3 = [maxlon_box, maxlat_box]; pt4 = [minlon_box, maxlat_box];
    pts = [pt1(1), pt2(1), pt3(1), pt4(1); pt1(2), pt2(2), pt3(2), pt4(2)];
    linpts = [pts'; pts(:, 1)'];
    m_line(linpts(:, 1), linpts(:, 2), 'color', 'm', 'linewidth', 2, 'linestyle', '-') 
    
    %Add mooring locations to map
    m_scatter(moorings.all(:, 2), moorings.all(:, 1), 250, 'k', 'p', 'filled')
    if splot == 1; m_scatter(moorings.all(:, 2), moorings.all(:, 1), 200, 'w', 'p', 'filled'); end
end

drawnow; pause(0.1)

if saveFigs 
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    saveas(gcf, [saveDir, saveName, '.fig'])
    print([saveDir, saveName],'-dpng')
end

clearvars -except rootPath ibcao AMSR2 profiles wvdata metData

