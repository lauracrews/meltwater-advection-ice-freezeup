%For a transect repeatedly sampled by the Seagliders, calculates latitude-bin-averaged
%observed mixed layer properties and concurrent modeled mixed layer
%properties. Calculates mixed layer heat budget for each transect. 

close all; 
clearvars -except rootPath AMSR2 profiles wvdata metData

saveFigs = true;

defineSODAconstants

%Heat content integration parameters
integrationDepth = 50; %Integrate to this depth [m] when calculating pycnocline heat content
densChange = 0.25; %Define base of mixed layer as density change [kg/m^3]

%Set up latitude bins in which to average observations
nbins = 8;
binInc = 0.8/nbins;
edgeLats = [73.6:binInc:74.4];
midLats = (edgeLats(1:end-1) + edgeLats(2:end))./2;

clim_salt = [25.5, 27]; %Color axis limits for plotting salinity on maps

%Set up map projection for the entire SODA-A to SODA-B area - used for
%making maps, but NOT for identifying profiles
minlon = -149; maxlon = -144; minlat = 72.85; maxlat = 74.55;%75.55; %
m_proj('lambert', 'lon', [minlon maxlon], 'lat', [minlat maxlat]);

%Geographic bounds for profiles on the southern Seaglider transect
minlat = 73 + 25/60; maxlat = 74.33; minlon = -148; maxlon = -146;

%Identify profiles in the geographic range of the Seaglider transect
inRegionMask = zeros(size(profiles.times));    
inRegionMask(profiles.lats >= minlat & profiles.lats <= maxlat ...
    & profiles.lons >= minlon & profiles.lons <= maxlon) = 1;

%%
plotCount = 1;

%Used to identify for export for the archive by write_dataForArchive.m
allProfs = [];
allTransects = [];

% Iterate through the Seaglider transects
for transect = 1:4

    disp(['Creating mixed layer heat budget for transect ', num2str(transect), ' of 4'])
    
    switch transect
        case 1 %First Seaglider transect           
            startTime = datenum('sept 26 2018') + .5;
            endTime = datenum('sept 30 2018');
        case 2 %Second Seaglider transect
            startTime = datenum('sept 30 2018');
            endTime = datenum('oct 6 2018');
        case 3 %Third Seaglider transect
            startTime = datenum('oct 6 2018');
            endTime = datenum('oct 11 2018')+ .2;
        case 4 %Fourth Seaglider transect
            startTime = datenum('oct 11 2018')+.2;
            endTime = datenum('oct 15 2018'); 
        case 5 %Western uCTD transect
            startTime = datenum('oct 11 2018')+.2;
            endTime = datenum('oct 15 2018'); 
    end
    
    %Indentify profiles in the time range
    inTimeMask = zeros(size(profiles.times));
    inTimeMask(profiles.times >= startTime & profiles.times < endTime) = 1;

    %For final time period, differentiate between Seaglider and uCTD profiles 
    platformMask = ones(size(profiles.times));   
    if transect == 4
        platformMask(strcmp(profiles.dataset, 'uCTD')) = 0;
    elseif transect == 5
        platformMask(~strcmp(profiles.dataset, 'uCTD')) = 0;
    end
    
    %Identify the profiles in this geographic region in this time period
    profNums = find(inTimeMask == 1 & inRegionMask == 1 & platformMask == 1 & profiles.qualFlag == 1); %inRegionMask was made outside of this for-loop
    
    allProfs = [allProfs, profNums];
    allTransects = [allTransects, transect .* ones(size(profNums))];
%%  
    %Bin the profiles by latitude, return the bin number of each profile
    [~, ~, binNums] = histcounts(profiles.lats(profNums), edgeLats);
    
    %% Calculate entrainment timeseries and F_atm timeseries for the initial set of profiles
    %This is used later to calculate how much entrainment and F_atm occured between
    %the initial observation and a later observation. This produces bin-averaged timeseries
    if transect == 1
           
        timeInc = 1/24; %Time interval in days at which to calculate running total of entrainment and Fatm. Speeds up calculation to not calculate entrainment at every 10-minute model timestep     
        binMeanTimeseriesTimes = [startTime:timeInc:datenum('oct 15 2018')]; %Times at which to calculate entrainment
        binMeanFatmTimeseries = nan.*ones(length(binMeanTimeseriesTimes), length(midLats)); %One column for each latitude bin
        binMeanEntrainmentTimeseries = nan.*ones(length(binMeanTimeseriesTimes), length(midLats)); %One column for each latitude bin
                
        %Empty arrays to hold the timeseries of Fatm and pycnocline heat
        %flux at each timestep in binMeanTimeseriesTimes
        FatmTimeseries_curProfiles = nan.*ones([length(binMeanTimeseriesTimes), length(profNums)]); 
        entrainmentTimeseries_curProfiles = nan.*ones([length(binMeanTimeseriesTimes), length(profNums)]); 
        
        for i = 1:length(profNums)
            profNum = profNums(i);
            
            % Mixed layer heat content of the current observed profile
            [mlHeatContent_obs, mlDepth_obs] = calculateHeatContent(profiles.CT(:, profNum), profiles.SA(:, profNum), profiles.sigthes(:, profNum), profiles.z, profiles.lats(profNum),...
                3, densChange); %Integrate to mixed layer depth

            % Heat content between the surface and integrationDepth of the current observed profile
            [upperHeatContent_obs, ~] = calculateHeatContent(profiles.CT(:, profNum), profiles.SA(:, profNum), profiles.sigthes(:, profNum), profiles.z, profiles.lats(profNum),...
                1, integrationDepth); %Integrate to a constant depth

            %Pycnocline heat content for the current observed profile
            initPycHeatContent = upperHeatContent_obs - mlHeatContent_obs;           

            %Load PWP results associated with the current observed profile
            [modelOutput, forc] = extract_pwpDataFromArchive(profNum);

            %Iterate through the time steps in the model output for the
            %current observed profile. Match with the times we want to
            %calculate entrainment
            [~, firstTimeInd] = min(abs(binMeanTimeseriesTimes - modelOutput.time(1))); %Begin filling in the timeseries at this binMeanTimeseriesTimes index - the closets binMeanTimeseriesTimes timestep to the observation time
            for curTimeInd = firstTimeInd:length(binMeanTimeseriesTimes) %curTimeInd is the time index in binMeanTimeseriesTimes
                [~, modelInd] = min(abs(binMeanTimeseriesTimes(curTimeInd) - modelOutput.time)); %Find the model output closes to this binMeanTimeseriesTimes timestep

                %Mixed layer heat content of the model output at the entrainment timestep  
                [mlHeatContent_pwp, ~] = calculateHeatContent(modelOutput.temp(:, modelInd), modelOutput.sal(:, modelInd), gsw_rho(modelOutput.sal(:, modelInd), modelOutput.temp(:,modelInd), 0), modelOutput.z, profiles.lats(profNum),...
                    3, densChange);

                %Upper ocean heate content of the model at the entrainment timestep
                [upperHeatContent_pwp, ~] = calculateHeatContent(modelOutput.temp(:, modelInd), modelOutput.sal(:, modelInd), gsw_rho(modelOutput.sal(:, modelInd), modelOutput.temp(:, modelInd), 0), modelOutput.z, profiles.lats(profNum),...
                    1, integrationDepth);

                modelPycHeatContent = upperHeatContent_pwp - mlHeatContent_pwp;

                %Will later average by diving entrainment (summed for
                %all profiles) by entrainmentCount (number of profiles
                %with entrainment data at this entraimentTimes timestep
                entrainmentTimeseries_curProfiles(curTimeInd, i) = initPycHeatContent - modelPycHeatContent;

                % Cumulative heat flux to the atmosphere from ERA5 data used to force this PWP run
                dt = forc.time(2)-forc.time(1); %Time increment between model timesteps, in days
%                 Fnet = forc.rhf + forc.thf; %Net surface heat flux used to force the model;
                Fnet = forc.heatFlux; %Net surface heat flux used to force the model;
                FatmTimeseries_curProfiles(curTimeInd, i) = (dt * 24*60*60) * trapz(-1 .* Fnet(1:modelInd)) / 10^6; %convert days to seconds, convert J to MJ        
            end
        end
        %At the end of these loops, have now created timeseries for each profile on the transect of
        %cumulative Fatm and cumulative heat flux from the pycnocline at each of the standard
        %timesteps in binMeanTimeseriesTimes. Now need to average these timeseries by latitude bin
        
        %Iterate through each latitude bin, average the entrainment and Fatm in that bin
        for binNum = 1:length(midLats)
            inds = find(binNums == binNum); %Indices of profiles in the current latitude bin          
            binMeanFatmTimeseries(:, binNum) = nanmean(FatmTimeseries_curProfiles(:, inds), 2);
            binMeanEntrainmentTimeseries(:, binNum) = nanmean(entrainmentTimeseries_curProfiles(:, inds), 2);
        end             
    end
   
    %% Calculate the entrainment between the times of observations on the current transect and 
    %the time of modelled freeze up for the PWP runs for those profiles
    
    %This differs from the previous entrainment calculation along the first
    %transect in that here we do not need to create a timeseries of
    %cumlative entrainment, just calculate the total entrainment bettween
    %the observation and the freeze up time
    observedMLheatContent = nan .* ones(size(profNums));
    observedUpperHeatContent = nan .* ones(size(profNums));
    observedMLdepth = nan .* ones(size(profNums));
    observedMLtemp = nan .* ones(size(profNums));
    observedMLsalt = nan .* ones(size(profNums));
    pwpMLheatContent_freezeTime = nan .* ones(size(profNums));
    pwpUpperHeatContent_freezeTime = nan .* ones(size(profNums));
    pwpMLdepth_freezeTime = nan .* ones(size(profNums));
    Fatm_obsToFreeze = nan .* ones(size(profNums));
    for i = 1:length(profNums)
        profNum = profNums(i);             
        
        %Mixed layer heat content of observed profile
        [mlHeatContent_obs, mlDepth_obs] = calculateHeatContent(profiles.CT(:, profNum), profiles.SA(:, profNum), profiles.sigthes(:, profNum), profiles.z, profiles.lats(profNum),...
            3, densChange);

        observedMLheatContent(i) = mlHeatContent_obs;
        observedMLdepth(i) = mlDepth_obs;
        
        %Average observed tempereature and salinity in the mixed layer
        [~, mlDepthInd] = min(abs(profiles.z - mlDepth_obs));
        observedMLtemp(i) = nanmean(profiles.CT(1:mlDepthInd, profNum));
        observedMLsalt(i) = nanmean(profiles.SA(1:mlDepthInd, profNum));      
        
        %Upper ocean heat conent of the observed profile
        [upperHeatContent_obs, ~] = calculateHeatContent(profiles.CT(:, profNum), profiles.SA(:, profNum), profiles.sigthes(:, profNum), profiles.z, profiles.lats(profNum),...
            1, integrationDepth);

        observedUpperHeatContent(i) = upperHeatContent_obs;
        
        %Calculate the above paramenters for the PWP profile at time of freeze-up
        densProfPWP = gsw_rho(profiles.pwpFreezeSAprof(:, profNum), profiles.pwpFreezeCTprof(:, profNum), 0);
        [mlHeatContent_pwp, mlDepth_pwp] = calculateHeatContent(profiles.pwpFreezeCTprof(:, profNum), profiles.pwpFreezeSAprof(:, profNum), densProfPWP, profiles.pwpz, profiles.lats(profNum),...
            3, densChange);
        
        pwpMLheatContent_freezeTime(i) = mlHeatContent_pwp;
        pwpMLdepth_freezeTime(i) = mlDepth_pwp;
        
        [upperHeatContent_pwp, ~] = calculateHeatContent(profiles.pwpFreezeCTprof(:, profNum), profiles.pwpFreezeSAprof(:, profNum), densProfPWP, profiles.pwpz, profiles.lats(profNum),...
            1, integrationDepth);

        pwpUpperHeatContent_freezeTime(i) = upperHeatContent_pwp;
               
        %%
         %Load the full PWP model run for the current observation
        [modelOutput, forc] = extract_pwpDataFromArchive(profNum);

        % Cumulative heat flux to the atmosphere from ERA5 data used to force this PWP run
        [~, endInd] = min(abs(forc.time - profiles.pwpFreezeTime(profNum)));
        dt = forc.time(2)-forc.time(1); %Time increment between model timesteps, in days
%         Fnet = forc.rhf + forc.thf; %Net surface heat flux used to force the model;
        Fnet = forc.heatFlux;
        Fatm_obsToFreeze(i) = (dt * 24*60*60) * trapz(-1 .* Fnet(1:endInd)) / 10^6; %convert days to seconds, convert J to MJ        
        
    end
    observedPycHeatContent = observedUpperHeatContent - observedMLheatContent; %Array of pycnocline contents at the times of the initial observations
    pwpPycHeatContent_freezeTime = pwpUpperHeatContent_freezeTime - pwpMLheatContent_freezeTime; %Array of PWP pycnocline contents at the times of the freeze up
    entrainment_obsToFreeze = observedPycHeatContent - pwpPycHeatContent_freezeTime; %Array of changes in pycnocline heat content = entrainment for each profile   
    depthChanges_obsToFreeeze = observedMLdepth - pwpMLdepth_freezeTime; %Array of changes in mixed layer depth from first observation until freeze up
    
    %% For later transects, calculate the amount of net surface heat flux and entrainment that happened
    %between the beginning of the study period and the currenet observation
    meanFatm_initToObs = zeros(size(profNums));
    meanEntrainment_initToObs = zeros(size(profNums));
    if transect > 1 %No need to do this for group 1, since can just use the model results from initial observation until freeze up
        for i = 1:length(profNums) %Iterate through all profiles on this transect
            profNum = profNums(i); %Current profile number        
            binNum = binNums(i); %Identify the bin the current profile belongs to              
            
            %Find the closest time in the entainment/Fatm timeseries created
            %from the original profiles to the time of the current profile.
            [~, endInd] = min(abs(binMeanTimeseriesTimes - profiles.times(profNum)));

            %The heat flux between the initial observations in the bin and the current observation
            meanFatm_initToObs(i) = binMeanFatmTimeseries(endInd, binNum);
            
            %The entrainment timeseries fills in the gap between the start
            %of the study period and the current profile
            meanEntrainment_initToObs(i) = binMeanEntrainmentTimeseries(endInd, binNum); %How much cumulative enrtrainment has happened in this bin since the initial time period
        end
    end     
            
    %% Add the entrainment from the initial time to the observation to the entrainment from the observation to freeze up
    entrainment_initToFreeze = meanEntrainment_initToObs + entrainment_obsToFreeze;
    Fatm_initToFreeze = meanFatm_initToObs + Fatm_obsToFreeze;
    
    %% Bin and average data    
    meanTimes = nan .* ones(size(midLats));
    meanFreezeUpTime = nan .* ones(size(midLats));
    stdFreezeUpTime = nan .* ones(size(midLats));

    meanObservedMLheat = nan .* ones(size(midLats));
    stdObservedMLheat = nan .* ones(size(midLats));
    
    meanObservedMLdepth = nan .* ones(size(midLats));
    stdObservedMLdepth = nan .* ones(size(midLats));

    meanObservedMLtemp = nan .* ones(size(midLats));
    stdObservedMLtemp = nan .* ones(size(midLats));
    meanObservedMLsalt = nan .* ones(size(midLats));
    stdObservedMLsalt = nan .* ones(size(midLats));
    
    meanPwpMLheatContent_freezeTime = nan .* ones(size(midLats));
    stdPwpMLheatContent_freezeTime = nan .* ones(size(midLats));
    
    meanEntrainment_initToFreeze = nan .* ones(size(midLats));
    stdEntrainment_initToFreeze = nan .* ones(size(midLats));
    
    meanFatm_initToFreeze = nan .* ones(size(midLats));
    stdFatm_initToFreeze = nan .* ones(size(midLats));
    for binNum = 1:length(midLats)
        curProfNums = profNums(binNums == binNum);
        
        meanTimes(binNum) = mean(profiles.times(curProfNums));
        meanFreezeUpTime(binNum) = mean(profiles.pwpFreezeTime(curProfNums));
        stdFreezeUpTime(binNum) = std(profiles.pwpFreezeTime(curProfNums));
        
        meanObservedMLheat(binNum) = mean(observedMLheatContent(binNums == binNum));
        stdObservedMLheat(binNum) = std(observedMLheatContent(binNums == binNum));
        meanObservedMLdepth(binNum) = mean(observedMLdepth(binNums == binNum));
        stdObservedMLdepth(binNum) = std(observedMLdepth(binNums == binNum));
        
        meanObservedMLtemp(binNum) = mean(observedMLtemp(binNums == binNum));
        stdObservedMLtemp(binNum) = std(observedMLtemp(binNums == binNum));
        meanObservedMLsalt(binNum) = mean(observedMLsalt(binNums == binNum));
        stdObservedMLsalt(binNum) = std(observedMLsalt(binNums == binNum));
        
        meanPwpMLheatContent_freezeTime(binNum) = mean(pwpMLheatContent_freezeTime(binNums == binNum));
        stdPwpMLheatContent_freezeTime(binNum) = std(pwpMLheatContent_freezeTime(binNums == binNum));
         
        meanEntrainment_initToFreeze(binNum) = mean(entrainment_initToFreeze(binNums == binNum));
        stdEntrainment_initToFreeze(binNum) = std(entrainment_initToFreeze(binNums == binNum));
        
        meanFatm_initToFreeze(binNum) = mean(Fatm_initToFreeze(binNums == binNum));
        stdFatm_initToFreeze(binNum) = std(Fatm_initToFreeze(binNums == binNum));
    end
        
   %% Find forward-running PWP model time from the model initialized in the earliest transect. 
   % Compare this concurrent model prediction with the observations taken on this later tranect
   %This is used to compare observed and modelled mixed layer
   %temperature, salinity, and depth. It is also used to calculate
   %advective heat loss
    meanPWPheat_obsTime = nan .* ones(size(midLats));
    stdPWPheat_obsTime = nan .* ones(size(midLats));
    meanPWPdepth_obsTime = nan .* ones(size(midLats));
    stdPWPdepth_obsTime = nan .* ones(size(midLats));  
    meanPWPmlTemp_obsTime = nan .* ones(size(midLats));
    stdPWPmlTemp_obsTime = nan .* ones(size(midLats));
    meanPWPmlSalt_obsTime = nan .* ones(size(midLats));
    stdPWPmlSalt_obsTime = nan .* ones(size(midLats));
   if transect > 1 %Because we want to compare later observations to the model results initiated on the first transect
              
       %Iterate through latitude bins, identify profiles from the current transect that are in the current latitude bin
       for binNum = 1:length(midLats)
            curProfNums = profNums(binNums == binNum);
            
           % Use this to compare everything against the earliest transect -
           % otherwise will compare against the most recent transect
            previousProfNums = group1ProfNums;

            %Bin the earlier profiles by latitude, identify which of the
            %previous profiles were in the current latitude bin
            [~, ~, prevBinNums] = histcounts(profiles.lats(previousProfNums), edgeLats);
            prevProfNums_curBin = previousProfNums(prevBinNums == binNum);
    
            %No current profiles in this bin at this timestep, but want an approximate time at which to plot
            %the mixed layer model results from this bin started at an
            %earlier timestep. Extrapolate from average time of profiles
            %taken in the adjacent bins
            if isempty(curProfNums) 
                meanTimes = fillmissing(meanTimes, 'nearest');
            end
                 
            %Empty arrays to hold PWP model predictions 
            pwpMLheatContent_obsTime = nan .* ones(size(prevProfNums_curBin));
            pwpMLdepth_obsTime = nan .* ones(size(prevProfNums_curBin)); 
            pwpMLtemp_obsTime = nan .* ones(size(prevProfNums_curBin)); 
            pwpMLsalt_obsTime = nan .* ones(size(prevProfNums_curBin)); 
            for i = 1:length(prevProfNums_curBin)
                
                %Load the full PWP model run for the observation from the
                %first transect
                [modelOutput, forc] = extract_pwpDataFromArchive(prevProfNums_curBin(i));                    

                %Closest pwp time step to the single observation
                [~, curTimeInd] = min(abs(modelOutput.time - meanTimes(binNum)));
                tempProf = modelOutput.temp(:, curTimeInd);
                saltProf = modelOutput.sal(:, curTimeInd);                    
                presProf = gsw_p_from_z(-modelOutput.z, modelOutput.config.lat); 
                densProfPWP = gsw_rho(saltProf, tempProf, 0);
                
                %Mixed layer heat content predicted by PWP at the time of
                %the observation on the current transect
                [mlHeatContent_pwp_obsTime, mlDepth_pwp_obsTime] = calculateHeatContent(tempProf, saltProf, densProfPWP, modelOutput.z, modelOutput.config.lat,...
                    3, densChange);
                
                pwpMLheatContent_obsTime(i) = mlHeatContent_pwp_obsTime;
                pwpMLdepth_obsTime(i) = mlDepth_pwp_obsTime;
                pwpMLtemp_obsTime(i) = tempProf(1);
                pwpMLsalt_obsTime(i) = saltProf(1);
                
            end
            %Average within latitude bins 
            meanPWPheat_obsTime(binNum) = nanmean(pwpMLheatContent_obsTime);
            stdPWPheat_obsTime(binNum) = std(pwpMLheatContent_obsTime);
            
            meanPWPdepth_obsTime(binNum) = nanmean(pwpMLdepth_obsTime);
            stdPWPdepth_obsTime(binNum) = std(pwpMLdepth_obsTime);
            
            meanPWPmlTemp_obsTime(binNum) = nanmean(pwpMLtemp_obsTime);
            stdPWPmlTemp_obsTime(binNum) = std(pwpMLtemp_obsTime);

            meanPWPmlSalt_obsTime(binNum) = nanmean(pwpMLsalt_obsTime);
            stdPWPmlSalt_obsTime(binNum) = std(pwpMLsalt_obsTime);
    
        end
    end
   
    advection = meanPWPheat_obsTime - meanObservedMLheat;
        
    %Save data from first transect for later comparisons.
    if transect == 1
        group1ProfNums = profNums;
        group1startTimes = meanTimes;
        group1mlHeats = meanObservedMLheat;
    end
    
    %Average time elapsed since the earliest observations in this bin,
    %Average heat fluxes - this is listed in the text of the document
    dt = (meanFreezeUpTime - group1startTimes) .* 24 .* 60 .* 60; %Convert days to seconds
    meanFatm_flux = (meanFatm_initToFreeze .* 10^6) ./ dt; %Convert MJ to J
    meanEntrainment_flux = (meanEntrainment_initToFreeze .* 10^6) ./ dt;

%     disp(['Transect ', num2str(transect), ' elapsed time until freeze up (by bin) = ', num2str(meanFreezeUpTime - group1startTimes), ' (days)'])
%     disp(['Transect ', num2str(transect), ' mean heat flux to atmosphere (by bin) = ', num2str(meanFatm_flux), ' (W/m^2)'])
%     disp(['Transect ', num2str(transect), ' mean heat flux from pycnocline (by bin) = ', num2str(meanEntrainment_flux), ' (W/m^2)'])
    
    %Print out average for whole transect, referenced in document text
%     disp(['Transect ', num2str(transect), ' mean F_atm = ', num2str(nanmean(Fatm_initToFreeze))])
%     disp(['Transect ', num2str(transect), ' mean F_pyc = ', num2str(nanmean(entrainment_initToFreeze))])
%     disp(['Transect ', num2str(transect), ' mean freeze up = ', datestr(nanmean(profiles.pwpFreezeTime(profNums)))])
%     if transect == 3; disp(['Transect ', num2str(transect), ' mean advection = ', num2str(nanmean(advection))]); end

    %% Begin plotting   
    figure(1); set(gcf, 'pos', [103 100 1519 827], 'color', 'w')
    
    %Map of 0-5 m salinity at profile locations
    for figNum = 1:2
        figure(figNum)
        switch figNum
            case 1
                subplot(4, 5, transect*5-4)
            case 2
                subplot(4, 3, transect*3-2)
        end
        h2 = m_scatter(profiles.lons(profNums), profiles.lats(profNums), 60, nanmean(profiles.salts(1:5, profNums)), 'filled');
        hold on
        
        %Do title now since it depends on the profiles numbers of the current transect
        title(['Observations ' datestr(min(profiles.times(profNums)), 'mmm dd'), ' to ', datestr(max(profiles.times(profNums)), 'mmm dd')], 'fontsize',12) 

    end
    
    figure(1)
    %Plot heat contents
    subplot(4, 5, transect*5 - 3); hold on
    
    %Concurrent PWP predicted heat content
    l(2) = bar(midLats, meanPWPheat_obsTime, 'FaceColor', 'none', 'LineWidth', 1);
    errorbar(midLats, meanPWPheat_obsTime, stdPWPheat_obsTime, 'LineStyle', 'none', 'LineWidth', 1, 'color', 'k')

    %Observed heat content
    l(1) = bar(midLats, meanObservedMLheat, .5, 'facealpha', 0.5);
    errorbar(midLats, meanObservedMLheat, stdObservedMLheat, 'LineStyle', 'none', 'LineWidth', 2, 'color', blue)
    ylabel('ML heat content (MJ/m^2)', 'fontsize', 12)
    ylim([0, 180]); xlim([edgeLats(1) - 0.05, edgeLats(end) + 0.05])
    xlabel('Latitude', 'fontsize', 12)
         
    lgd = legend(l(1:2), 'Observed', '1D prediction');
    set(lgd, 'fontsize', 12, 'location', 'north')    
    
    % Plot of observed and modelled mixed layer temperature
    subplot(4, 5, transect*5 - 2); hold on
    m(1) = scatter(midLats, meanObservedMLtemp, 40, blue, 'filled');
    errorbar(midLats, meanObservedMLtemp, stdObservedMLtemp,'LineStyle', 'none', 'LineWidth', 1, 'color', blue)
    m(2) = scatter(midLats, meanPWPmlTemp_obsTime, 40, 'k', 'filled');
    errorbar(midLats, meanPWPmlTemp_obsTime, stdPWPmlTemp_obsTime,'LineStyle', 'none', 'LineWidth', 1, 'color', 'k')
    ylabel(['Mixed layer temperature (', sprintf(char(176)), 'C)'], 'fontsize', 12)
    ylim([-1.5, 0])
    lgd = legend(m(1:2), 'Observed', '1D Prediction');
    set(lgd, 'fontsize', 12, 'location', 'north') 
%     disp(['Transect ', num2str(transect), ' observed mixed layer temperature = ', num2str(nanmean(observedMLtemp))])
%     disp(['Transect ', num2str(transect), ' modeled-observed mixed layer temperature = ', num2str(nanmean(meanPWPmlTemp_obsTime-meanObservedMLtemp))])
    
    % Plot of observed and modelled mixed layer salinity
    subplot(4, 5, transect*5 - 1); hold on
    scatter(midLats, meanObservedMLsalt, 40, blue, 'filled')
    errorbar(midLats, meanObservedMLsalt, stdObservedMLsalt,'LineStyle', 'none', 'LineWidth', 1, 'color', blue)
    scatter(midLats, meanPWPmlSalt_obsTime, 40, 'k', 'filled')
    errorbar(midLats, meanPWPmlSalt_obsTime, stdPWPmlSalt_obsTime,'LineStyle', 'none', 'LineWidth', 1, 'color', 'k')
    ylabel('Mixed layer salinity (g/kg)', 'fontsize', 12)
    ylim([25, 26.8])
%     disp(['Transect ', num2str(transect), ' observed mixed layer temperature = ', num2str(nanmean(observedMLsalt))])

    % Plot of observed and modelled mixed layer depth
    subplot(4, 5, transect*5); hold on
    scatter(midLats, meanObservedMLdepth, 40, blue, 'filled')
    errorbar(midLats, meanObservedMLdepth, stdObservedMLdepth,'LineStyle', 'none', 'LineWidth', 1, 'color', blue)
    scatter(midLats, meanPWPdepth_obsTime, 40, 'k', 'filled')
    errorbar(midLats, meanPWPdepth_obsTime, stdPWPdepth_obsTime,'LineStyle', 'none', 'LineWidth', 1, 'color', 'k')
    ylabel('Mixed layer depth (m)', 'fontsize', 12)
    ylim([5, 35])
%     disp(['Transect ', num2str(transect), ' modeled mixed layer depth = ', num2str(nanmean(meanPWPdepth_obsTime))])
%     disp(['Transect ', num2str(transect), ' observed mixed layer depth = ', num2str(nanmean(observedMLdepth))])
%     disp(['Transect ', num2str(transect), ' modeled-observed mixed layer depth = ', num2str(nanmean(meanPWPdepth_obsTime-meanObservedMLdepth))])

    figure(2);
    
    %Plot of modelled freeze up times
    subplot(4, 3, transect*3-1); hold on
    scatter(midLats, meanFreezeUpTime, 40, 'k', 'filled')    
    errorbar(midLats, meanFreezeUpTime, stdFreezeUpTime,'LineStyle', 'none', 'LineWidth', 1, 'color', 'k')
    ylim([datenum('oct 4, 2018'), datenum('oct 29, 2018')])
    ylabel('Modeled freeze up time', 'fontsize', 12)
    datetick('y',  'mmm-dd', 'keeplimits', 'keepticks')
    
    %Plot of heat budget components
    subplot(4, 3, transect*3);
    totalHeatLoss = [meanFatm_initToFreeze; -meanEntrainment_initToFreeze; advection]';
        
    br = bar(midLats, totalHeatLoss);    
    clrs = [blue; [0 0 0]; yellow];
    set(br, {'facecolor'}, {clrs(1,:),clrs(2,:), clrs(3,:)}.')
    if transect > 1
        lgd = legend('Atmosphere', 'Pycnocline', 'Advection');
        set(lgd,  'fontsize', 12)
    end
    ylabel('Heat out of mixed layer (MJ/m^2)', 'fontsize', 12)
    ylim([-50, 150])
        
    plotCount=plotCount+1;

end

%% Format and save figures
for figNum = 1:2 
    switch figNum
        case 1 %Figure showing mixed layer properties
            numCols = 5;
            ps = [103 82 1310 873];
            saveDir = [rootPath, 'figures/fig8/'];
            saveName = 'observed_modeled_mixedLayerProperties';
        case 2 %Figure showing modeled freeze up date and heat budget components
            numCols = 3;
            ps = [445 82 800 871];
            saveDir = [rootPath, 'figures/fig9/'];
            saveName = 'mixedLayerHeatBudgets';
    end
    figure(figNum)
    set(gcf, 'color', 'w', 'pos', ps)

    for splot = 1:(4*numCols)
       subplot(4, numCols, splot); 
       if mod(splot, numCols) == 1 %Format maps
           cmocean('haline'); caxis(clim_salt) %Constant color axis everywhere
            m_grid('xtick', [-148:2:-144], 'fontsize', 12)   
            m_scatter(moorings([1, 3:4], 2), moorings([1, 3:4], 1), 250, 'k', 'p', 'filled')
        %     m_scatter(moorings(:, 2), moorings(:, 1), 200, 'w', 'p', 'filled')
       else
            xlabel(['Latitude (', sprintf(char(176)), 'N)'], 'fontsize', 12)
            set(gca, 'XTick', edgeLats, 'XTickLabelRotation', 90)
            xlim([edgeLats(1) - 0.05, edgeLats(end) + 0.05])
            grid on; box on
       end
    end
    
    if saveFigs 
        if ~exist(saveDir, 'dir'); mkdir(saveDir); end
        print([saveDir, saveName],'-dpng')
        saveas(gcf, [saveDir, saveName, '.fig'])
    end
end
       
        