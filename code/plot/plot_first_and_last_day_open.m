%Plot map of melt-out and freeze-up date from AMSR2
%Adapted from code from Luc Rainville

close all;
clearvars -except rootPath AMSR2 profiles wvdata metData

%Load bathymetry
load('ibcao.mat')

dayThreshold = 3; %Number of days in a fow that that a cell must be ice covered or ice free
[moorings, colors] = defineSODAconstants; %Mooring locations and standard colors
% load('Beaufort_AMSR2_2018.mat') %Load AMSR2 data
sz = size(AMSR2.lon);

saveFigs = true;
saveDir = [rootPath ,'figures/fig1/'];
saveName = 'ice_formation_melt';

%% Calculate days of melt out and freeze up

% first day when the ice concentration is >0.15 for more than DT days.
meltOut = zeros(size(AMSR2.lon)); %length(bremen.lat), length(bremen.lon))*NaN;

% last day when the ice concentration is >0.15 for more than DT days.
freezeUp = zeros(size(AMSR2.lon)); %zeros(length(bremen.lat), length(bremen.lon))*NaN;

% mininum ice (was it open at all during the summer?)
minIce = zeros(size(AMSR2.lon)); %zeros(length(bremen.lat), length(bremen.lon));

%Calculate first and last day (or detect if always ice covered) for each
%lat/lon grid point
for k1=1:sz(1)
    for k2=1:sz(2)
        ii = conv2P(double(squeeze(AMSR2.SIC(k1,k2,:))<0.15), ones(dayThreshold,1)/dayThreshold); % this is one only if it was true DT/2 days before and after
        if ~isempty(find(ii==1, 1))
            meltOut(k1,k2)=AMSR2.mattime(find(ii==1,1));
            freezeUp(k1,k2)=AMSR2.mattime(find(ii==1,1,'last'));
        else 
            minIce(k1,k2)=1; % completely ice covered the whole time. 
        end
    end
end

%%
lonlimits = [-170 -120];
latlimits = [68 80];
m_proj('Lambert','lon',lonlimits,'lat',latlimits)
year = 2018;

%Limits for colormaps
colorIncrement_meltout = 7; %step size for colormap, in days
colorLims_meltOut = [datenum(year,7,15) datenum(year,9,15)]; %Earliest and latest day for color axis
colorIncrement_freezeUp = 3; %step size for colormap, in days
colorLims_freezeUp = [datenum(year,10,1) datenum(year,10,30)];

%Begin plotting melt out timing
meltOut_plot = meltOut;
meltOut_plot(meltOut < colorLims_meltOut(1)) = colorLims_meltOut(1);
meltOut_plot(minIce == 1) = nan;

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
freezeUp_plot(freezeUp<colorLims_freezeUp(1))=colorLims_freezeUp(1);
freezeUp_plot(minIce == 1) = nan;

h2=subplot_m(1, 2, 2, 0.05);
m_pcolor(AMSR2.lon, AMSR2.lat, freezeUp_plot); shading flat;
caxis([colorLims_freezeUp(1)-colorIncrement_freezeUp/2 colorLims_freezeUp(2)+colorIncrement_freezeUp/2]);
title({'Date of ice formation',num2str(year)},'fontsize',20)

%Set up freeze up colormap
c_tick = colorLims_freezeUp(1):colorIncrement_freezeUp:colorLims_freezeUp(2)+colorIncrement_freezeUp/2;
colorLims_freezeUp=[colorLims_freezeUp(1) c_tick(end)];
cmap = cbrewer('seq','Greys', length(c_tick));

%Make colorbar for freeze up time
c1 = colorbar_out(0.5,0, 0.02);
clear c_tick_label; count = 1;
for k=1:2:length(c_tick) %Create vector of colorbar labels, labelling every other color increment
  c_tick_label{count} = datestr(c_tick(k),'mmm dd');
  count = count + 1;
end
pp = get(c1, 'position');
set(c1,'ytick', c_tick(1:2:end), 'YTickLabel',c_tick_label, 'position',[pp(1)-pp(3),pp(2)+pp(4)/3,pp(3:4)]);
colormap(h2,cmap)

%Additions to both plots
for splot = 1:2
    h2=subplot_m(1, 2, splot, 0.05);
    hold on
    m_contour(ibcao.LON, ibcao.LAT,ibcao.d,[-1000 -1000],'color',[1 1 1]*0.9,'linewidth',1); %Plot bathymetry
    m_contour(AMSR2.lon, AMSR2.lat, minIce,[1 1]*0.15,'color',[0 0 1],'linewidth',1)
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
    print([saveDir, saveName],'-dpng')
    saveas(gcf, [saveDir, saveName, '.fig'])
end
