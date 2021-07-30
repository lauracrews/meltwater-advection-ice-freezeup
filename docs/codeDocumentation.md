## Data extraction and processing

#### AMSR2 daily sea ice concentrations
`batchAMSR2.m` - from Samuel Brenner. Uses ftp to download daily sea ice concentration data from the University of Bremen archive 

`readAMSR2.m` - from Samuel Brenner. Reads daily .hdf files of sea ice concentration into the Matlab data structure `AMSR2_2018.mat`
##

#### CryoSat-2 dynamic ocean topography
`extract_DOT.m` - Reads CyroSat-2 .asc files into Matlab and extracts lat/lon and dynamic ocean topography (DOT) - some string manipulation required. Calculates time of each data point, which are referenced to the satellite's equator crossing. Saves `allDOT.mat` data structure combining all data and optionally saves .mat data structures corresponding to each individual .asc file

## Make figures

#### Figure 1
`plot_first_and_last_day_open.m` - Makes Figure 1, maps of melt out and freeze up date. Adapted from code from Luc Rainville. Requires the `AMSR2_2018.mat` data structure, which was created by running `batchAMSR2.m` and `readAMSR2.m`. Uses the function `conv2P.m` from Luc Rainville

#### Figure 2
`plot_surfaceTemperature_byTime.m` - Plot the surface ocean temperature measurements from *Healy*, Seagliders, uCTD divided by time periods

`plot_surfaceSalinity_byTime.m` - Plot the surface salinity measurements from *Healy*, Seagliders, uCTD divided by time periods

#### Figure 3
`plot_observations_dailyMODIS.m` - Compares in situ sea surface temperature to MODIS sea surface temperature averaged in a box surrounding the observation. Iterates through daily MODIS images to find the closest MODIS SST in time to each observation, up to a threshold number of days away from the observation day. Options to plot data from all platforms on one figure or to make separate panels for the *Healy* underway data, Seaglider and uCTD data, and Wave Glider data.

`plot_observations_SMOSsalinity.m` - Makes plot of SMOS sea surface salinity along with sea surface salinity from in situ observations. Iterates through weekly SMOS images and finds in situ observations taken at the time of those images. Options to plot data from all platforms on one figure or to make separate panels for *Healy* underway data, Seaglider and uCTD data, and Wave Glider data. 

#### Figure 4
`plot_remnantIceSAR.m` - Plots the sea surface salinity from Seaglider and Wave Glider in the remnant ice area overlaid on concurrent SAR imagery.

`plot_frontCharacteristics.m` - Plots crossing of meltwater front by Wave Glider and Seaglider. Calculates and plots horizontal buoyancy gradient, vertical buoyancy frequency, geostrophic velocity, Richardson number, and Rossby number. The subplots for the density and the horizontal buoyancy gradient along the Waveglider track are used in the document. Other calculations are discussed in the text. 

`plot_averageProfiles.m` - Identifies profiles taken within a specified geographic region and time period. Plots those profiles vs. depth as well as the average profile. Used in Figure 4 to plot the profiles taken in the remnant ice area.

#### Figure 5
`plot_section.m` - Plots temperature and salinity sections from Seagliders and uCTD as well as maps showing the locations of the profiles.

#### Figure 6
`plot_profileTS.m` - Makes temperature-salinity plots for profiles within a given geographic
and time range. Also makes a map showing where the profiles were taken.

#### Figure 7
`plot_DOT_geostrophic.m` - Calculates geostrophic currents from the CryoSat-2 dynamic ocean topography (extracted to Matlab using `extract_DOT.m`). Plots DOT along with geostrophic current vectors. 

`plot_frontTracking.m` - Creates and plots simulated trajectories. Makes a vector time series of Ekman velocities. Optionally animates the trajectories 

#### Figure 8
`plot_heatBudgetTransects.m` - For a transect repeatedly sampled by the Seagliders, calculates latitude-bin-averaged observed mixed layer properties and concurrent modeled mixed layer properties. Calculates mixed layer heat budget for each transect

#### Figure 9
Also made with `plot_heatBudgetTransects.m` described above. The heat budget calculations are outlined in the flowchart shown below
![](https://github.com/lauracrews/meltwaterAdvection/blob/main/docs/Heat%20Budget%20Methods.png)
*This flowchart outlines the method for creating heat budgets for each Seaglider transect between the time of the first Seaglider transect and the modeled freeze up time* 

#### Figure 10
`plot_pwpResultsComparison.m` - Plots initial temperature and salinity profiles as well as the resulting PWP profiles at time of freezing for specified regions and time periods. Uses a specified example profile for each region/time


## Miscellaneous calculations
`calculate_averageProperties_atFreezeUp.m` - calculates the average heat content in the upper 40 m at modeled freeze up for profiles inside and outside the meltwater - these values are used in the discussion section of the paper.

`calculateCompassAngle.m` - Takes in zonal and meridional velocity vector components and calculates the compass angle of that velocity (i.e. north = 0°, east = 90°). Used in `plot_frontTracking.m`

`calculateHeatContent.m` - Depth-integrates the heat content relative to freezing of a T-S profile. Uses salinity data to calculate freezing point. If surface data is missing, fills with the shallowest measurement. Then if additional data is missing fills with linear interpolation. 

`defineSODAconstants.m` - Various useful values defined in one place: standard colors used throughout, locations of the moorings within the study area
