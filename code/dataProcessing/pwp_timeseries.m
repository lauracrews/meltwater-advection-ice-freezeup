%Initializes and runs the PWP model code. Intakes data structure of
%atmospheric data with fields for surface stress, heat fluxes, precip/evap.
%Initializes PWP with temperature and  salinity profiles that are the  
%average of the startProfNums indicies in the the matrix of input profiles 

%Note: Laura's code for using ice growth for freshwater forcing is in pwp_timeseries_draft.m
%It has been deleted here since the questions adressed using PWP are
%related to the time it takes to reach freezing, not the evolution of
%the mixed layer after ice begins to form
function [forc, config, init, modelOutput, atmFluxes] = pwp_timeseries(mldThreshold, pwpMaxDepth, atmFluxes, profiles, startProfNums, useFreshwaterForcing, integrationType, integrationLimit)
        
    warning('off','MATLAB:interp1:NaNinY')
    clear model init config forc    
    
    %Configuration
    config.lat = nanmean(profiles.lats(startProfNums)); %latitude
    config.Kz = 1e-6; %Scalar vertical diffusion
    config.Km = 1e-6; %Momentum vertical diffusion
    config.r= 0;%;1.2662e-6; %damping term
    config.tintv= 6; %interval at which the resutls are plotted
    config.dz = 1; %depth increment (meters)
    config.zmax= pwpMaxDepth; %the depth to run
    config.dt= 10*60; %time-step increment (convert minutes to seconds)
    config.I0= [0.6200 0.3800]; %[0.6200 0.3800]: Eq. 6 of PWP86
    config.lambda= [0.6000 20]; %[0.6000 20], longwave (0.6 m) and shortwave (20 m) extinction coefficients
    config.rho= 1026; %density of seawater 
    config.cp= cp_t_exact(27, -1, 0); %specific heat of water (J/kgC)
    config.BRiCrit= 0.65; %critical bulk richardson number (0.65)
    config.GRiCrit= 0.3; %critical gradient richardson number (0.25)
    config.plotflag= 0;
    config.writeout= 1; %PWP code crashes if you turn this off... just saves the file, you can delete it later
    config.runname = 'pwpOut';% [profiles.dataset{startProfNums(1)}, '_profile', num2str(startProfNums(1)), 'to_', profiles.dataset{startProfNums(end)}, '_profile', num2str(startProfNums(end))];
    
    % Forcing
    startTime = nanmean(profiles.times(startProfNums));
    startInd = find([atmFluxes.time] >= startTime, 1);
    atmTime = [atmFluxes(startInd:end).time]'; %Time of atmospheric forcing
    forc.time = atmTime(1):(config.dt/(60*60*24)):atmTime(end); %Convert time step in seconds to fractional days   
    
    %Added by Laura to make it easy to vary the MLD definition used in pwp_MixedLayerDepth.m 
    forc.mldThreshold = mldThreshold;   
    
    %Interpolate freshwater forcing to the timestep used in PWP
    if useFreshwaterForcing
        netFWFlux = ([atmFluxes(startInd:end).mtpr]' + [atmFluxes(startInd:end).mer]').* 1e-3; %Note that evaporation is already negative. Convert mm/sec to m/sec 
        forc.FWFlux = interp1(atmTime, netFWFlux , forc.time); 
    else
        forc.FWFlux = zeros(size(forc.time));
    end
     
    %Interpolate heat and surface stress forcing to the timestep used in PWP
    forc.rhf = interp1(atmTime, [atmFluxes(startInd:end).msnswrf]', forc.time); %Net shortwave 
    netSurfaceHeat = [atmFluxes(startInd:end).msnlwrf]' + [atmFluxes(startInd:end).msshf]' + [atmFluxes(startInd:end).mslhf]';
    forc.thf = interp1(atmTime, netSurfaceHeat, forc.time);%Heat added to surface        
    
    forc.taux = interp1(atmTime, [atmFluxes(startInd:end).metss]', forc.time);
    forc.tauy = interp1(atmTime, [atmFluxes(startInd:end).mntss]', forc.time);

    %Note: Laura's code for adding a simulated storm to the surface stress
    %forcing is in pwp_timeseries_draft.m

    %Initial CTD profile
    init.z  = (config.dz/2:config.dz:config.zmax-config.dz/2)';  % Vector of depth values
    
    %If running PWP on the average of multiple profiles, this will average
    %them together. If running on one profile, no effect
    tempProf = nanmean(profiles.CT(:, startProfNums), 2);
    saltProf = nanmean(profiles.SA(:, startProfNums), 2);
    
    %Interpolate to the PWP depth grid
    tempProf = interp1(profiles.z, tempProf, init.z);
    saltProf = interp1(profiles.z, saltProf, init.z);
    
    %If surface data is missing, fill with the shallowest available mreasurement
    %This was done in loadProfiles.m by it's necessary to repeat  
    %because NaN is introduced to the top cell when interpolating to the PWP depth grid
    firstTemp = find(~isnan(tempProf), 1, 'first');
    tempProf(1:firstTemp) = tempProf(firstTemp);    
    firstSalt = find(~isnan(saltProf), 1, 'first');
    saltProf(1:firstSalt) = saltProf(firstSalt);

    %If measurements don't go down to 50 m, interpolate down
    lastTemp = find(~isnan(tempProf), 1, 'last');
    tempProf(lastTemp:end) = tempProf(lastTemp);    
    lastSalt = find(~isnan(saltProf), 1, 'last');
    saltProf(lastSalt:end) = saltProf(lastSalt);
       
    %Profiles to initialize PWP
    init.S = saltProf;
    init.T = tempProf;
    init.UV = zeros(length(init.z),2);
    init.startProfNums = startProfNums;

    %Run the model
    modelOutput = pwp_Run(config,init,forc);
    
    %Iterate through the PWP output and calculate the heat content of each profile
    for pwpProfNum = 1:length(modelOutput.time)
        tempProfPWP = modelOutput.temp(:, pwpProfNum);
        saltProfPWP = modelOutput.sal(:, pwpProfNum);
        densProfPWP = gsw_rho(saltProfPWP, tempProfPWP, 0);
        modelOutput.pden(:, pwpProfNum) = densProfPWP;

        [heatContent, integrationDepth] = calculateHeatContent(tempProfPWP, saltProfPWP, densProfPWP, modelOutput.z, config.lat, integrationType, integrationLimit); 
        modelOutput.heatContent(pwpProfNum) = heatContent;
        modelOutput.integrationDepth(pwpProfNum) = integrationDepth;
    end
end      
