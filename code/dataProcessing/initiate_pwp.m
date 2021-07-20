%This script identifies observed temperature/salinity profiles in a region 
%and time period and initializes and runs the pwp model on each of those
%profiles, after averaging ERA5 near the profile to force model. Saves the 
%pwp model output at each timestep for each profile, as well as adding key
%pwp results (i.e. profiles at freezing and freeze up time) to the
%profiles.mat data structure. 
%(PWP runs multiple times - slow!)
close all; clear

integrationDepthType = 1; %Integrate heat content to a constant depth
integrationDepth = 40; %depth to integrate to
mldThreshold = 0.01; %Change in density relative to surface [kg/m^3]. Used in pwp_MixedLayerDepth to determine the depth range over which to distribute momentum
pwpMaxDepth = 50;
useFreshwaterForcing = true;
if useFreshwaterForcing %Where model output for individual profiles will be saved
    savedir = '/Users/lcrews/Documents/MATLAB/SODA/finalFigures/pwpResults/withFreshwater';
else
    savedir = '/Users/lcrews/Documents/MATLAB/SODA/finalFigures/pwpResults/noFreshwater';
end

%Use data from between these dates
startTime = datenum('sept 21 2018'); %Earliest profiles to use
endTime_profiles = datenum('oct 31 2018'); %Latest profiles to use
endTime_model = datenum('nov 7 2018'); %Run PWP models until this time

%Geographic bounds to identifying in situ profiles - same as boundaries
%in plot_surfaceSalinity_byTime.m
minlon = -149; maxlon = -143; minlat = 72.75; maxlat = 75.55; 
pt1 = [minlon, minlat]; pt2 = [maxlon, minlat]; pt3 = [maxlon, maxlat]; pt4 = [minlon, maxlat];
pts = [pt1(1), pt2(1), pt3(1), pt4(1); pt1(2), pt2(2), pt3(2), pt4(2)];

%Identify profiles in this region and time period. profilesInRegion_timeseries 
%will call a function to integrate their heat content
profiles = profilesInRegion_timeseries(pts, startTime, endTime_profiles, integrationDepthType, integrationDepth);  

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
    [forc, config, init, modelOutput] = pwp_timeseries(mldThreshold, pwpMaxDepth, atmFluxes, profiles, i, useFreshwaterForcing, useConstantDepth, integrationDepth);                
    cd(savedir); save(['pwp_initProfile', num2str(i), '_Output.mat'], 'forc', 'config', 'init', 'modelOutput');
    
    %Find which PWP profile had the upper ocean reach the freezing point
    pwpTimestep = 1; freezeInd = [];
    while pwpTimestep <= length(modelOutput.time)       
        freezingPtProf = gsw_CT_freezing(modelOutput.sal(:, pwpTimestep), gsw_p_from_z(-modelOutput.z, profiles.lats(profNum)));       
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

%Save the entire profiles datastructure, used later for plotting
cd('/Users/lcrews/Documents/MATLAB/SODA/finalFigures/pwpResults')
if useFreshwaterForcing
    save('profiles_withFreshwater.mat', 'profiles') 
    save('profiles_boundaries_withFreshwater.mat', 'pts')
else
    save('profiles_noFreshwater.mat', 'profiles')
    save('profiles_boundaries_noFreshwater.mat', 'pts')
end