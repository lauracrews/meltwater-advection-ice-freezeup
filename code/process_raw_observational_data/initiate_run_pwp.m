%This script identifies observed temperature/salinity profiles in a region 
% and time period and initializes and runs the pwp model on each of those
% profiles, after averaging ERA5 near the profile to force model. Saves the 
% pwp model output at each timestep for each profile, as well as adding key
% pwp results (i.e. profiles at freezing and freeze up time) to the
% profiles.mat data structure. 
%(PWP runs multiple times - slow!)
function profiles = initiate_run_pwp(profiles, startTime, endTime_model, integrationType, integrationDepth, useFreshwaterForcing, pwpMaxDepth)

    
    mldThreshold = 0.01; %Change in density relative to surface [kg/m^3]. Used in pwp_MixedLayerDepth to determine the depth range over which to distribute momentum
    if useFreshwaterForcing %Where model output for individual profiles will be saved
        savedir = '/Users/lcrews/Documents/MATLAB/meltwaterAdvection/data/pwpResults/withFreshwater/';
    else
        savedir =  '/Users/lcrews/Documents/MATLAB/meltwaterAdvection/data/pwpResults/noFreshwater/';
    end
    
    %Initialize fields in the "profiles" data structure to hold results from
    %the PWP model ("profiles" will just hold PWP results at time of freeze-up,
    %while the PWP data for the entire model run will be saved separately into
    %its own .mat file
    profiles.pwpFreezeTime = nan .* ones(size(profiles.lats));
    profiles.pwpFreezeHeatContent = nan .* ones(size(profiles.lats));
    profiles.pwpFreezeCTprof = nan .* ones(pwpMaxDepth, length(profiles.lats));
    profiles.pwpFreezeSAprof = nan .* ones(pwpMaxDepth, length(profiles.lats));
    profiles.pwpFreezePdenProf = nan .* ones(pwpMaxDepth, length(profiles.lats));
    
    %% Initialize and run the PWP model for each observed profile in the region and time period
    for i = 1:length(profiles.lats)
        
        disp(['Processing profile ', num2str(i), ' of ', num2str(length(profiles.lats))])
        
        %Calculate atmosphehric forcing in a box surrounding the current profile
        minlon = profiles.lons(i) - 0.5; maxlon = profiles.lons(i) + 0.5;  %At 75 N, 0.15 degrees latitude ~= 16 km, 0.5 degrees longitude ~= 14.5 km
        minlat = profiles.lats(i) - 0.15; maxlat = profiles.lats(i) + 0.15; 
        pt1 = [minlon, minlat]; pt2 = [maxlon, minlat]; pt3 = [maxlon, maxlat]; pt4 = [minlon, maxlat];
        pts = [pt1(1), pt2(1), pt3(1), pt4(1); pt1(2), pt2(2), pt3(2), pt4(2)];
        atmFluxes = makeERA5timeseries(pts, startTime, endTime_model, {'msnlwrf', 'msnswrf', 'mslhf', 'msshf', 'metss', 'mntss', 'mer', 'mtpr'});
    
        %Run PWP on the current profile (profile i) and save PWP output
        [forc, config, init, modelOutput] = pwp_timeseries(mldThreshold, pwpMaxDepth, atmFluxes, profiles, i, useFreshwaterForcing, integrationType, integrationDepth);                
        save([savedir, 'pwp_initProfile', num2str(i), '_Output.mat'], 'forc', 'config', 'init', 'modelOutput');
        disp(['Saving all pwp reults for profile ', num2str(i), ' as ', savedir, 'pwp_initProfile', num2str(i), '_Output.mat'])

        %Find which PWP profile had the upper ocean reach the freezing point
        pwpTimestep = 1; freezeInd = [];
        while pwpTimestep <= length(modelOutput.time)       
            freezingPtProf = gsw_CT_freezing(modelOutput.sal(:, pwpTimestep), gsw_p_from_z(-modelOutput.z, profiles.lats(i)));       
            if modelOutput.temp(3, pwpTimestep) <= freezingPtProf(3)
                freezeInd = pwpTimestep;
                break
            end
            pwpTimestep = pwpTimestep + 1;
        end
    
        if ~isempty(freezeInd) %Otherwise defaults to the nan it was initialized with
            profiles.pwpFreezeTime(i) = modelOutput.time(freezeInd);
            profiles.pwpFreezeHeatContent(i) = modelOutput.heatContent(freezeInd);
            profiles.pwpFreezeCTprof(:, i) = modelOutput.temp(:, freezeInd);
            profiles.pwpFreezeSAprof(:, i) = modelOutput.sal(:, freezeInd);
            profiles.pwpFreezePdenProf(:, i) = modelOutput.pden(:, freezeInd);
        end
        
    end
    profiles.pwpz = modelOutput.z;
    delete('pwpOut_Output.mat') %This is created automatically by the PWP script, we don't want it
    
end