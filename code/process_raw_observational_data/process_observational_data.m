% This a master script for loadding, processing, and saving datasets that
% were uploaded to UW Researchworks at https://digital.lib.washington.edu/researchworks/items/e5410a95-d046-4458-989d-e87c327d8e59
% 
% It completes various tasks with the observational data:

% -Corrects the temperature and salinity data from Healy's underway system
%
% -Loads Waveglider data, writes to .nc file (minimal processin)
%
% -Creates combined data structure of all uCTD an Seaglider profiles in the
% region of interest. Calculates heat content of these. Checks quality
% 
% -nitiates the PWP model with observed temperature andd salinity
% profiles. Runs the PWP model with ERA5 atmospheric forcing near the
% observation used to initialize the model

%% Where to save the processed data
archiveddata_savedir =  '/Users/lcrews/Documents/MATLAB/meltwaterAdvection/dataForArchive/';

%% Make a .csv file with the data collected by Healy's seawater intake system
%while underway
healyUnderway = load_correct_underwayData;
writetable(struct2table(healyUnderway), [archiveddata_savedir, 'healyUnderway_test.csv'])

%% Make a .nc file of the Wave Glider data. Code adapted from SWIFT2NC.m by
%Jim Thomson
waveglider = loadWaveglider;
filename = 'waveglider.nc';
write_glider_ncfile([archiveddata_savedir, filename], waveglider);

%% Make a .nc file of the Seaglider and uCTD profile data

% Specify time and geographic bounds for uCTD or Seaglider profiles

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
integrationDepthType = 1; %Integrate heat content to a constant depth
integrationDepth = 40; %depth to which to integrate heat
profiles = profilesInRegion_timeseries(pts, startTime, endTime_profiles, integrationDepthType, integrationDepth);  

%% Run PWP model on observed profiles

% Running the model is slow - we don't want to rerun later. 
% Run once and then save the profiles datastructure including model results
run_model = false;
datadir = '/Users/lcrews/Documents/MATLAB/meltwaterAdvection/data/';
useFreshwaterForcing = true; 

if run_model

    % Note! The "profiles" data structure will just hold PWP results at time of freeze-up,
    % while the full PWP results for the each model run will be saved separately into
    % its own .mat file in ./meltwaterAdvection/data/pwpResult/

    pwpMaxDepth = 50;
    profiles = initiate_run_pwp(profiles, startTime, endTime_model, integrationDepthType, integrationDepth, useFreshwaterForcing, pwpMaxDepth);
    
    if useFreshwaterForcing
        save([datadir, 'profiles_pwpresults_withFreshwater.mat'], 'profiles') 
    else
        save([datadir, 'profiles_pwpresults_noFreshwater.mat'], 'profiles')
    end
else %Model already run

    if useFreshwaterForcing
        load([datadir, 'profiles_pwpresults_withFreshwater.mat'])
    else
        load([datadir, 'profiles_pwpresults_noFreshwater.mat'])
    end
end

%Save profiles including PWP model results to a .nc file
filename = 'profiles.nc';
write_profiles_ncfile([archiveddata_savedir, filename], profiles);

%% Write separate .nc files for the PWP model results for each model run
% used in the heat budget

modelresults_savedir = [archiveddata_savedir, 'pwpResults/'];
if ~exist(modelresults_savedir, 'dir')
    mkdir(modelresults_savedir)
end

%Contains .mat files listing the profile ID number of each observational
%profile on each of the four heat budget transects
heatbudget_datadir = '/Users/lcrews/Documents/MATLAB/meltwaterAdvection/data/heat_budget_transects/';
files = dir([heatbudget_datadir, '*.mat']);
all_heat_budget_profiles = [];
for filenum = 1:length(files)
    load(files(filenum).name, 'profNums')
    all_heat_budget_profiles = [all_heat_budget_profiles; profNums(:)];
end

datadir = '/Users/lcrews/Documents/MATLAB/meltwaterAdvection/data/';
if useFreshwaterForcing
    pwp_results_dir = [datadir, 'pwpResults/withFreshwater/'];
else
    pwp_results_dir = [datadir, 'pwpResults/noFreshwater/']
end

for i = 1:length(all_heat_budget_profiles)
    profNum = all_heat_budget_profiles(i);   
     
    %Load PWP data for this profile
    load([pwp_results_dir, 'pwp_initProfile', num2str(profNum), '_Output.mat'], 'modelOutput', 'forc')
    
    %Begin making netCDF of PWP data for this profile
    filename = ['pwpResults_profile', num2str(profNum), '.nc'];
    
    write_pwp_model_ncfile([modelresults_savedir, filename], modelOutput, forc, profNum, profiles);
end


%% Functions to write .nc files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = write_glider_ncfile(filename, waveglider)
    ncid=netcdf.create(filename,'CLOBBER');
    t_dim = netcdf.defDim(ncid,'time', length(waveglider.times));
    z_dim = netcdf.defDim(ncid,'z', size(waveglider.z, 2));
    names=fieldnames(waveglider);
    
    %Create the variables
    for i = 1:length(names) 
        if i < 5
            eval(strcat(names{i}, '_id = netcdf.defVar(ncid, ''', names{i}, ''', ''NC_DOUBLE'', [t_dim]);'));
        else
            eval(strcat(names{i}, '_id = netcdf.defVar(ncid, ''', names{i}, ''', ''NC_DOUBLE'', [t_dim z_dim]);'));  
        end
    end
    netcdf.endDef(ncid);
    
    %Fill in values for the variables
    for i = 1:length(names)
        eval(strcat('netcdf.putVar(ncid, ', names{i}, '_id, ', 'waveglider.', names{i}, ');'));  
    end  
    netcdf.close(ncid)
    
    %Add units and descriptions
    ncwriteatt(filename,'vehicle','units','none')
    ncwriteatt(filename,'vehicle','long_name','Wave Glider ID number')
    
    ncwriteatt(filename,'lons','units','degrees E')
    ncwriteatt(filename,'lons','long_name','Longitude')
    
    ncwriteatt(filename,'lats','units','degrees N')
    ncwriteatt(filename,'lats','long_name','Latitude')
    
    ncwriteatt(filename,'times','units','Days from Jan 0, 0000')
    ncwriteatt(filename,'times','long_name','Time of measurement')
    
    ncwriteatt(filename,'z','units','m')
    ncwriteatt(filename,'z','long_name','Depth of CTD below ocean surface (positive downward)')
    
    ncwriteatt(filename,'temps','units','degrees C')
    ncwriteatt(filename,'temps','long_name','Water temperature')
    
    ncwriteatt(filename,'salts','units','psu')
    ncwriteatt(filename,'salts','long_name','Salinity')
    
    ncwriteatt(filename,'CT','units','degrees C')
    ncwriteatt(filename,'CT','long_name','Conservative temperature')
    
    ncwriteatt(filename,'SA','units','g/kg')
    ncwriteatt(filename,'SA','long_name','Absolute salinity')
    
    ncwriteatt(filename,'sigthes','units','kg m^-3')
    ncwriteatt(filename,'sigthes','long_name','Potential density referenced to the surface')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write .nc file for profiles data 
function [] = write_profiles_ncfile(filename, profiles)
    
    vehicleNames = nan .* ones(size(profiles.dataset));
    for j = 1:length(profiles.dataset)
        curName = profiles.dataset{j};
        if strcmp(curName, 'uCTD')
            vehicleNames(j) = 1;
        else
            vehicleNames(j) = str2double(curName(3:end));
        end
    end
    profiles.vehicle = vehicleNames;
    profiles.pwpProfNum = 1:length(profiles.lats);
    
    %Rename some fields
    profiles = renameStructField(profiles, 'pwpFreezeCTprof', 'pwpCT');
    profiles = renameStructField(profiles, 'pwpFreezeSAprof', 'pwpSA');
    profiles = renameStructField(profiles, 'pwpFreezePdenProf', 'pwpSigthes');
    
    %Remove derived variables that we don't want to save
    profiles = rmfield(profiles, {'dataset', 'integrationDepths', 'heatContents', 'pwpFreezeHeatContent'});
    
    %Reorder remaining fields to be in the same order as other data sets
    P = [15 16 2 3 4 1 5 7 6 8 9 17 10 14 11 12 13]; 
    profiles = orderfields(profiles, P);
    %%
    %Begin making netCDF
    ncid=netcdf.create(filename,'CLOBBER');
    prof_dim = netcdf.defDim(ncid,'profileNumber', length(profiles.times));
    z_dim = netcdf.defDim(ncid,'z', size(profiles.z, 1));
    pwpz_dim = netcdf.defDim(ncid,'pwp_z', size(profiles.pwpz, 1));
    
    names=fieldnames(profiles);
    
    %Create the variables
    for i = 1:length(names)
        
        %Define the dimensions of each variable
        if i <= 5 || i == 12 || i == 13
            dimText = '[prof_dim]';
        elseif i == 6
            dimText = '[z_dim]';
        elseif i >= 7 && i <= 11
            dimText = '[z_dim prof_dim]';
        elseif i == 14
            dimText = '[pwpz_dim]';
        elseif i >= 15
            dimText = '[pwpz_dim prof_dim]';
        end
        
        eval(strcat(names{i}, '_id = netcdf.defVar(ncid, ''', names{i}, ''', ''NC_DOUBLE'',', dimText, ');'));
    end
    netcdf.endDef(ncid);
    
    %Fill in values for the variables
    for i = 1:length(names)
        eval(strcat('netcdf.putVar(ncid, ', names{i}, '_id, ', 'profiles.', names{i}, ');'));
    end  
    netcdf.close(ncid)
    
    %Add units and descriptions
    ncwriteatt(filename,'vehicle','units','none')
    ncwriteatt(filename,'vehicle','long_name','Seaglider ID# or, if equal to 1, indicates uCTD')
    
    ncwriteatt(filename,'lons','units','degrees E')
    ncwriteatt(filename,'lons','long_name','Longitude')
    
    ncwriteatt(filename,'lats','units','degrees N')
    ncwriteatt(filename,'lats','long_name','Latitude')
    
    ncwriteatt(filename,'times','units','Days from Jan 0, 0000')
    ncwriteatt(filename,'times','long_name','Time of measurement')
    
    ncwriteatt(filename,'z','units','m')
    ncwriteatt(filename,'z','long_name','Depth below ocean surface (positive downward)')
    
    ncwriteatt(filename,'temps','units','degrees C')
    ncwriteatt(filename,'temps','long_name','Water temperature')
    
    ncwriteatt(filename,'salts','units','psu')
    ncwriteatt(filename,'salts','long_name','Salinity')
    
    ncwriteatt(filename,'CT','units','degrees C')
    ncwriteatt(filename,'CT','long_name','Conservative temperature')
    
    ncwriteatt(filename,'SA','units','g/kg')
    ncwriteatt(filename,'SA','long_name','Absolute salinity')
    
    ncwriteatt(filename,'sigthes','units','kg m^-3')
    ncwriteatt(filename,'sigthes','long_name','Potential density referenced to the surface')
    
    ncwriteatt(filename,'pwpz','units','m')
    ncwriteatt(filename,'pwpz','long_name','Center depth of PWP model depth bin (positive downward)')
    
    ncwriteatt(filename,'pwpProfNum','units','None')
    ncwriteatt(filename,'pwpProfNum','long_name','ID# of this profile to match to PWP model output in pwpResults.nc')
    
    ncwriteatt(filename,'pwpFreezeTime','units','days from Jan 0, 0000')
    ncwriteatt(filename,'pwpFreezeTime','long_name','Time the 3rd depth bin (centered at 2.5 m) of the PWP model fell below freezing')
    
    ncwriteatt(filename,'pwpCT','units','degrees C')
    ncwriteatt(filename,'pwpCT','long_name','PWP modeled conservative temperature at pwpFreezeTime')
    
    ncwriteatt(filename,'pwpSA','units','g/kg')
    ncwriteatt(filename,'pwpSA','long_name','PWP modeled absolute salinity at pwpFreezeTime')
    
    ncwriteatt(filename,'pwpSigthes','units','kg m^-3')
    ncwriteatt(filename,'pwpSigthes','long_name','PWP modeled potential density at pwpFreezeTime')
    
    ncwriteatt(filename,'qualFlag','units','True/False')
    ncwriteatt(filename,'qualFlag','long_name','Flag = 1 if profile considered "good", flag = 0 if excluded from analysis')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
function [] = write_pwp_model_ncfile(filename, modelOutput, forc, profNum, profiles)
    ncid=netcdf.create(filename,'CLOBBER');
    
    %Define array dimensions
    t_dim = netcdf.defDim(ncid,'time', length(forc.time));
    z_dim = netcdf.defDim(ncid,'z', size(modelOutput.z, 1));
    
    %Create the variables
    z_id = netcdf.defVar(ncid, 'z', 'NC_DOUBLE', [z_dim]);
    time_id = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', [t_dim]);
    heat_id = netcdf.defVar(ncid, 'heatFlux', 'NC_DOUBLE', [t_dim]);
    temp_id = netcdf.defVar(ncid, 'CT', 'NC_DOUBLE', [z_dim t_dim]);
    salt_id = netcdf.defVar(ncid, 'SA', 'NC_DOUBLE', [z_dim t_dim]);
    sigthe_id = netcdf.defVar(ncid, 'sigthe', 'NC_DOUBLE', [z_dim t_dim]);

    netcdf.endDef(ncid);

    %Fill the variables
    netcdf.putVar(ncid, z_id, modelOutput.z);
    netcdf.putVar(ncid, time_id, forc.time);
    netcdf.putVar(ncid, heat_id, forc.rhf + forc.thf);
    netcdf.putVar(ncid, temp_id, modelOutput.temp);
    netcdf.putVar(ncid, salt_id, modelOutput.sal);
    netcdf.putVar(ncid, sigthe_id, modelOutput.pden);
        
    netcdf.close(ncid)

    %Define global attributes
    ncwriteatt(filename, '/', 'profileNumber', profNum);
    ncwriteatt(filename,'/', 'heatBudgetTransect', nan);%allTransects(i));
    ncwriteatt(filename,'/', 'lat_observedProfile', profiles.lats(profNum));
    ncwriteatt(filename,'/', 'lon_observedProfile', profiles.lons(profNum));
    ncwriteatt(filename,'/', 'time_observedProfile', datestr(profiles.times(profNum), 'mmm dd yyyy hh:MM:ss'));
        
    %Attributes of the variables 
    ncwriteatt(filename,'time','units','Days from Jan 0, 0000')
    ncwriteatt(filename,'time','long_name','Time of PWP model timestep')
    
    ncwriteatt(filename,'z','units','m')
    ncwriteatt(filename,'z','long_name','Center depth of PWP model depth bin (positive downward)')
   
    ncwriteatt(filename,'heatFlux','units', 'W m^-2')
    ncwriteatt(filename,'heatFlux','long_name', 'Net surface heat flux')
    
    ncwriteatt(filename,'CT','units', 'degrees C')
    ncwriteatt(filename,'CT','long_name', 'Conservative temperature of model output')
    
    ncwriteatt(filename,'SA','units', 'g/kg')
    ncwriteatt(filename,'SA','long_name', 'Absolute salinity of model output')
    
    ncwriteatt(filename,'sigthe','units','kg m^-3')
    ncwriteatt(filename,'sigthe','long_name', 'Potential density  of model output referenced to the surface')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%