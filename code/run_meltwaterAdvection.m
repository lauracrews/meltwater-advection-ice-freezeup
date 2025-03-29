% Reproduces the analysis and figures of "Direct observations of the role
% of lateral advection of sea ice meltwater in the onset of autumn freeze
% up" by Crews et al. (2022).

% All scripts by Laura Crews (lcrews@uw.edu) unless otherwise noted 
% July 2021.

% To get started, download the required data to the appropriate directories. 
% Instructions are provided in the "download data" section of 
% https://github.com/lauracrews/meltwaterAdvection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all

rootPath = [userpath, '/meltwaterAdvection/'];
saveFigs = true;

% If running observational data processing and pwp modeling:
% process_observational_data

% If data processing is complete (i.e., using archived data)
% Open the the Seaglider, underway CTD, Wave Glider, and USCGC Healy 
% underway data used in this study (download data from http://hdl.handle.net/1773/47135)
disp('Reading observational data into Matlab')
[metData, wvdata, profiles] = extract_dataFromArchive(rootPath);

%%

if ~exist([rootPath, 'data/AMSR2_2018.mat'], 'file')
    disp('Downloading and extracting AMSR2 sea ice concentration')
    batchAMSR2_1819(rootPath); %Download AMSR2 sea ice concentration from University of Bremen
    AMSR2 = readAMSR2; %Create and load .mat file from AMSR2 data
else 
    load AMSR2_2018.mat %Mat file has already been prepared. Just reload it
end


%Begin making figures
if ~exist([rootPath, 'figures/'], 'dir')
    mkdir([rootPath, 'figures/'])
end

%Load bathymetry for figures
load('ibcao.mat')

%%

disp('Making Figure 1 by running plot_first_and_last_day_open.m')
plot_first_and_last_day_open; 

%You need to download the MODIS-Terra SST data from the Physical Oceanography 
%Distributed Active Archive Center at NASA/JPL.
disp('Making Figure 2 (top panels) by running plot_surfaceTemperature_byTime.m')
plot_surfaceTemperature_byTime

disp('Making Figure 2 (bottom panels) by running plot_surfaceSalinity_byTime.m')
plot_surfaceSalinity_byTime

if ~exist([rootPath, 'data/modisComparison.mat'], 'file')
    disp('Matching observations to MODIS-SST data - this will take some time')
    match_observations_modis
end

%Note that Figure 3 won't be saved until after
%plot_observations_SMOSsalinity.m has been run
disp('Making Figure 3a by running plot_observations_dailyMODIS.m')
plot_observations_dailyMODIS
 
%You need to download the SMOS sea surface salinity data from SEANOE
disp('Making Figure 3b by running plot_observations_SMOSsalinity.m')
plot_observations_SMOSsalinity

%Note: need to run Figure 4 construction in this order to not have a bug
%that distorts the multiple map axes
disp('Making Figure 4e and 4f by running plot_averageProfiles.m')
plot_averageProfiles

disp('Making Figure 4a and 4d by running plot_remnantIceSAR.m')
plot_remnantIceSAR

disp('Making Figure 4b and 4c by running plot_frontCharacteristics.m')
plot_frontCharacteristics

disp('Making Figure 5 by running plot_section.m')
for transect = 1:2
    plot_section
end

disp('Making Figure 6 by running plot_profileTS.m')
plot_profileTS

%You need to download the dynamic ocean typography (DOT) data (available at https://github.com/lauracrews/meltwaterAdvection/blob/main/rawDOT.zip) 
%Extract the DOT data from .asc files to one large .mat file
if ~exist([rootPath, 'data/allDOT.mat'], 'file')
    extract_DOT
end

disp('Starting Figure 7 by running plot_DOT_geostrophic.m')
plot_DOT_geostrophic

disp('Finishing Figure 7 by running plot_frontTracking.m')
plot_frontTracking

%You need to download the PWP model results used to make the heat budget (download from http://hdl.handle.net/1773/47135)
disp('Making figures 8 and 9 by running plot_heatBudgetTransects.m')
plot_heatBudgetTransects

disp('Making Figure 10 by running plot_pwpResultsComparison.m')
plot_pwpResultsComparison
