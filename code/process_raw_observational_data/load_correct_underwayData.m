%Compares Healy underway temperature and salinity data to concurrent uCTD
%data. Applies a linear fit correction to the underway salinity data and a constant
%correction to the underway temperature data. 

%Reference: Alory et al., 2015, "The French contribution to the voluntary
%observing ships network of sea surface salinity", doi:10.1016/j.dsr.2015.08.005

function metData = load_correct_underwayData()
%%
    makePlots = false;
    timeThreshold = 30; %Minutes on either side of a uCTD measurement to include avereaged underway data for comparison. 1-hr total windows (so timeThreshold = 30) used in Alory et al. 2015
    timeThreshold = timeThreshold/(60*24); %Convert minutes to days

    %Load underway data from the Healy system, this data was previously read into a .mat file using load_HEALY1802_MET.m 
    load met.mat
    metTime = met.time; metLat = met.LA; metLon = met.LO; 
    metSeaTemp = met.TT; metSalt = met.SA;

    %Calculate conservative temperature and absolute salinity for the underway data
    [metSA, ~, ~] = gsw_SA_Sstar_from_SP(metSalt, gsw_p_from_z(-2.7.*ones(size(metLat)), metLat), metLon, metLat);
    metCT = gsw_CT_from_t(metSA, metSeaTemp, gsw_p_from_z(-2.7.*ones(size(metLat)), metLat));

    %Load all uCTD profiles (not just those in the study area). Eliminate
    %profiles with bad data (mostly due to salinity cell clogging)
    load uCTD_SODA_2018.mat; load goodProfiles.mat
    firstProf = 11; %Skip the profiles right at the shelf break
    goodProfiles = goodProfiles(firstProf:end);
    
    %"Near surface" uCTD temperature and salinity averaged in the upper 5 m
    uctdSalt = nanmean(uCTD.salinity(1:5, goodProfiles), 1); uctdSalt = uctdSalt(:);
    uctdTemp = nanmean(uCTD.temperature(1:5, goodProfiles), 1); uctdTemp = uctdTemp(:);
            
    %Calculate conservative temperature and absolute salinity for the uCTD data
    [uctdSA, ~, ~] = gsw_SA_Sstar_from_SP(uctdSalt', nanmean(uCTD.pressure(1:5, goodProfiles), 1), uCTD.lon(goodProfiles), uCTD.lat(goodProfiles));
    uctdCT = gsw_CT_from_t(uctdSA, uctdTemp', nanmean(uCTD.pressure(1:5, goodProfiles), 1));
    uctdSA = uctdSA(:); uctdCT = uctdCT(:);
    
    %Use elapesed time throughout the cruise instead of the very large
    %datetime format times for the linear regressions
    elapsedTimes_uCTD = uCTD.time(goodProfiles) - min(uCTD.time(goodProfiles));
    elapsedTimes_met = metTime - min(uCTD.time(goodProfiles));
    
    if makePlots %Set up map projection
        minlon = -156; maxlon = -136; minlat = 70; maxlat = 78;
        m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);
    end

    %Identify the underway measurements around each uCTD cast
    concurrentMetTemps = nan .* ones(size(goodProfiles));
    concurrentMetSalts = nan .* ones(size(goodProfiles));
    concurrentMetCTs = nan .* ones(size(goodProfiles));
    concurrentMetSAs = nan .* ones(size(goodProfiles));
    for i = 1:length(goodProfiles)  
        profNum = goodProfiles(i);
        
        %Time range of underway data to include
        [~, metStartInd] = min(abs(metTime - (uCTD.time(profNum)-timeThreshold)));
        [~, metEndInd] = min(abs(metTime - (uCTD.time(profNum)+timeThreshold)));

        %Median underway data 
        concurrentMetTemps(i) = nanmedian(metSeaTemp(metStartInd:metEndInd));
        concurrentMetSalts(i) = nanmedian(metSalt(metStartInd:metEndInd));
        
        concurrentMetCTs(i) = nanmedian(metCT(metStartInd:metEndInd));
        concurrentMetSAs(i) = nanmedian(metSA(metStartInd:metEndInd));

        %Optional code to plot the measurement locations and values  
        if makePlots
            m_scatter(uCTD.lon(profNum), uCTD.lat(profNum), 40, nanmean(uCTD.temperature(1:5, profNum)), 'filled');
            hold on
            m_scatter(metLon(metStartInd:metEndInd), metLat(metStartInd:metEndInd), 20, metCT(metStartInd:metEndInd), 'filled');
        end

    end
    
    if makePlots; m_grid; end

    %Eliminate a couple of outlier salinitiess
    notOut = find(abs(concurrentMetSAs - uctdSA) <= 1.5); %these are acceptable values

    %Linear fit of salinity differences over time
    mdlSalt = fitlm(elapsedTimes_uCTD(notOut), uctdSalt(notOut) - concurrentMetSalts(notOut));
    SaltCorrections = mdlSalt.Coefficients.Estimate(2) .* elapsedTimes_met + mdlSalt.Coefficients.Estimate(1);
    metSalt_corrected = metSalt + SaltCorrections;
    
    mdlSA = fitlm(elapsedTimes_uCTD(notOut), uctdSA(notOut) - concurrentMetSAs(notOut));
    SAcorrections = mdlSA.Coefficients.Estimate(2) .* elapsedTimes_met + mdlSA.Coefficients.Estimate(1);
    metSA_corrected = metSA + SAcorrections;
    
    %Linear fit of temperature differences over time
%     mdlTemp = fitlm(elapsedTimes_uCTD, uctdTemp - concurrentMetTemps);
%     TempCorrections = mdlTemp.Coefficients.Estimate(2) .* elapsedTimes_met + mdlTemp.Coefficients.Estimate(1);
%     metTemp_corrected = metSeaTemp + TempCorrections;
%     
%     mdlCT = fitlm(elapsedTimes_uCTD, uctdCT - concurrentMetCTs);
%     CTcorrections = mdlCT.Coefficients.Estimate(2) .* elapsedTimes_met + mdlCT.Coefficients.Estimate(1);
%     metCT_corrected = metCT + CTcorrections;
    
    %Constant temperature correction
    TempCorrection = median(uctdTemp - concurrentMetTemps);
    metTemp_corrected = metSeaTemp + TempCorrection .* ones(size(metSeaTemp));
    
    CTcorrection = median(uctdCT - concurrentMetCTs);
    metCT_corrected = metCT + CTcorrection .* ones(size(metSeaTemp));
    
    %Eliminate a bunch of NaNs at the start. This is cleaner for archiving
    %the data
    ind1 = find(~isnan(metSeaTemp), 1, 'first');
    ind2 = find(~isnan(metSeaTemp), 1, 'last');
    
    %Add all data to data structure to be returned by the function
    metData.lats = metLat(ind1:ind2); metData.lons = metLon(ind1:ind2); metData.times = metTime(ind1:ind2);
    metData.temps_uncorrected = metSeaTemp(ind1:ind2); metData.temps = metTemp_corrected(ind1:ind2);
    metData.CTs_uncorrected = metCT(ind1:ind2); metData.CTs = metCT_corrected(ind1:ind2); 
    metData.salts_uncorrected = metSalt(ind1:ind2); metData.salts = metSalt_corrected(ind1:ind2);
    metData.SAs_uncorrected = metSA(ind1:ind2); metData.SAs = metSA_corrected(ind1:ind2); 
    
    if makePlots
        figure; set(gcf, 'color', 'w')
        subplot(2, 1, 1)
        scatter(elapsedTimes_met, metData.SAs_uncorrected);
        hold on; grid on
        scatter(elapsedTimes_uCTD, uctdSA, 40, 'r', 'filled');
%         xlim([min(uCTD.time(goodProfiles)) - .5, max(uCTD.time(goodProfiles)) + .5])
%         datetick('x', 'keeplimits')
        lgd = legend('TSG S_A', 'uCTD S_A');
        set(lgd, 'fontsize', 12)
        ylim([22, 30])
        xlim([min(elapsedTimes_uCTD) - 0.5, max(elapsedTimes_uCTD) + .5])
        xlabel('Elapsed time since first uCTD cast (days)', 'fontsize', 12)


        subplot(2, 1, 2);
        scatter(elapsedTimes_uCTD(notOut), uctdSA(notOut) - concurrentMetSAs(notOut), 40, 'r', 'filled');
        hold on; grid on
        scatter(elapsedTimes_uCTD(out), uctdSA(out) - concurrentMetSAs(out), 40, 'k', 'filled');
        xlim([min(elapsedTimes_uCTD) - .5, max(elapsedTimes_uCTD) + .5])
        h = plot(mdlSA); delete(h([1, 3:4]));
        set(h(2), 'color' ,'k', 'linewidth', 1)
%         datetick('x', 'keeplimits')
        ylabel('uCTD S_A - TSG S_A', 'fontsize', 12)
        lgd = legend('Salinity difference', 'Excluded outlier', ['Linear Fit y = ', num2str(round(mdlSA.Coefficients.Estimate(2), 3)), 'x + ', num2str(round(mdlSA.Coefficients.Estimate(1), 2))]);
        set(lgd, 'fontsize', 12, 'location', 'southeast')
        title(''); %xlabel('');
        xlabel('Elapsed time since first uCTD cast (days)', 'fontsize', 12)
 
        %Temperature data
        figure; set(gcf, 'color', 'w')
        subplot(2, 1, 1)
        scatter(elapsedTimes_met, metData.CTs_uncorrected);
        hold on; grid on
        scatter(elapsedTimes_uCTD, uctdCT, 40, 'r', 'filled');
%         xlim([min(elapsedTimes_uCTD) - .5, max(elapsedTimes_uCTD) + .5])
%         datetick('x', 'keeplimits')
        lgd = legend('TSG CT', 'uCTD CT');
        set(lgd, 'fontsize', 12)
        ylim([-3, 3])
        xlim([min(elapsedTimes_uCTD) - 0.5, max(elapsedTimes_uCTD) + .5])
        xlabel('Elapsed time since first uCTD cast (days)', 'fontsize', 12)


        %Temeprature differences
        subplot(2, 1, 2); 
        scatter(elapsedTimes_uCTD, uctdCT - concurrentMetCTs, 40, 'r', 'filled');
        xlim([min(elapsedTimes_uCTD) - .5, max(elapsedTimes_uCTD) + .5])
        hold on; grid on
        h = plot(mdlCT); delete(h([1, 3:4]));
        set(h(2), 'color' ,'k', 'linewidth', 1)
        ylabel('uCTD CT - TSG CT', 'fontsize', 12)
        lgd = legend('Temperature difference', ['Linear Fit y = ', num2str(round(mdlCT.Coefficients.Estimate(2), 3)), 'x + ', num2str(round(mdlCT.Coefficients.Estimate(1), 2))]);
        set(lgd, 'fontsize', 12, 'location', 'southeast')
        title(''); %xlabel('');
%         datetick('x', 'keeplimits')
        xlabel('Elapsed time since first uCTD cast (days)', 'fontsize', 12)
    end



    if makePlots
        figure
        subplot(1, 2, 1)
        scatter(uctdCT, concurrentMetCTs, 30)%, elapsedTimes_uCTD)
        caxis([min(uCTD.time), max(uCTD.time)])
        cb = colorbar; datetick(cb, 'y', 'keeplimits');
        ylabel(cb, 'Measurement time', 'fontsize', 12)
        hold on; grid on
        mdl = fitlm(uctdCT, concurrentMetCTs);
        h = plot(mdl); delete(h([1, 3:4]));
        set(h(2), 'color' ,'r', 'linewidth', 1)
        limits = [-2, 2];
        xlim(limits), ylim(limits);
        x = limits(1):.1:limits(2);
        y = x; plot(x, y, 'color', 0.7 .* [1 1 1], 'LineStyle', '-')
        pbaspect([1 1 1])
        xlabel('uCTD conservative temperature avg 0-5 m', 'fontsize', 12);
        ylabel('Ship conservative temperature at 2.7 m', 'fontsize', 12)
        lgd = legend('Observed temperature', ['Linear Fit y = ', num2str(round(mdl.Coefficients.Estimate(2), 2)), 'x + ', num2str(round(mdl.Coefficients.Estimate(1), 2))]);
        set(lgd, 'fontsize', 12, 'location', 'southeast')
        title('');

        subplot(1, 2, 2)
        scatter(uctdSA, concurrentMetSAs, 30, elapsedTimes_uCTD)
        hold on
        % scatter(uctdSA(notOut), concurrentMetSAs(notOut), 50, 'k');
        scatter(uctdSA(out), concurrentMetSAs(out), 50, 'k', 'filled')
        caxis([min(uCTD.time), max(uCTD.time)])
        cb = colorbar; datetick(cb, 'y', 'keeplimits');
        ylabel(cb, 'Measurement time', 'fontsize', 12)
        hold on
        grid on
        mdl = fitlm(uctdSA(notOut), concurrentMetSAs(notOut));
        h = plot(mdl); delete(h([1, 3:4]));
        set(h(2), 'color' ,'r', 'linewidth', 1)
        limits = [24.5, 31];
        xlim(limits), ylim(limits);
        x = limits(1):.1:limits(2);
        y = x; plot(x, y, 'color', 0.7 .* [1 1 1], 'LineStyle', '-')
        pbaspect([1 1 1])
        xlabel('uCTD absolute salinity avg 0-5 m', 'fontsize', 12);
        ylabel('Ship absolute salinity at 2.7 m', 'fontsize', 12)

        lgd = legend('Observed salinity', 'Excluded outlier', ['Linear Fit y = ', num2str(round(mdl.Coefficients.Estimate(2), 2)), 'x + ', num2str(round(mdl.Coefficients.Estimate(1), 2))]);
        set(lgd, 'fontsize', 12, 'location', 'southeast')
        title('');
    end

end
 