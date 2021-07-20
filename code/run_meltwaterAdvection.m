%Reproduces the analysis and figures of "Meltwater advection hastens autumn
%freeze up" by Crews et al. (2021). 

%All scripts by Laura Crews (lcrews@uw.edu) unless otherwise noted 
%July 2021.

%First download the Seaglider, underway CTD, Wave Glider, and USCGC Healy 
%underway data used in this study to the directory ~/meltwaterAdvection/data/ 
%from the University of Washington ResearchWorks archive 
%at http://hdl.handle.net/1773/47135
disp('Reading observational data into Matlab')
extract_dataFromArchive

disp('Downloading and extracting AMSR2 sea ice concentration')
batchAMSR2; %Download AMSR2 sea ice concentration from University of Bremen
readAMSR2; %Create .mat file from AMSR2 data
load AMSR2_2018.mat

if ~exist([userpath, '/meltwaterAdvection/figures/'], 'dir')
    mkdir([userpath, '/meltwaterAdvection/figures/'])
end

%Need to make ibcao.mat available to contour bathymetry
disp('Making Figure 1 by running plot_first_and_last_day_open.m')
plot_first_and_last_day_open; 

%First download the MODIS-Terra SST data from the Physical Oceanography 
%Distributed Active Archive Center at NASA/JPL.
disp('Making Figure 2 (top panels) by running plot_surfaceTemperature_byTime.m')
plot_surfaceTemperature_byTime

disp('Making Figure 2 (bottom panels) by running plot_surfaceSalinity_byTime.m')
plot_surfaceSalinity_byTime

if ~exist([userpath, '/meltwaterAdvection/data/modisComparison.mat'], 'file')
    disp('Matching observations to MODIS-SST data - this will take some time')
    match_observations_modis
end

%Note that Figure 3 won't be saved until after
%plot_observations_SMOSsalinity.m has been run
disp('Making Figure 3a by running plot_observations_dailyMODIS.m')
plot_observations_dailyMODIS
 
%First download the SMOS sea surface salinity data from SEANOE
disp('Making Figure 3b by running plot_observations_SMOSsalinity.m')
plot_observations_SMOSsalinity

disp('Making Figure 4a and 4d by running plot_remnantIceSAR.m')
plot_remnantIceSAR

disp('Making Figure 4b and 4c by running plot_frontCharacteristics.m')
plot_frontCharacteristics

disp('Making Figure 4e and 4f by running plot_averageProfiles.m')
plot_averageProfiles

disp('Making Figure 5 by running plot_section.m')
for transect = 1:2
    plot_section
end

disp('Making Figure 6 by running plot_profileTS.m')
plot_profileTS

%First download the dynamic ocean typography (DOT) data. 
%Extract the DOT data from .asc files to one large .mat file
cd([userpath, '/meltwaterAdvection/data/'])
if ~exist('allDOT.mat', 'file')
    extract_DOT
end

disp('Starting Figure 7 by running plot_DOT_geostrophic.m')
plot_DOT_geostrophic

disp('Finishing Figure 7 by running plot_frontTracking.m')
plot_frontTracking

disp('Making figures 8 and 9 by running plot_heatBudgetTransects.m')
plot_heatBudgetTransects

disp('Making Figure 10 by running plot_pwpResultsComparison.m')
plot_pwpResultsComparison
