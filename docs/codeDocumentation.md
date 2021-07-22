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

## Miscellaneous calculations
`calculate_averageProperties_atFreezeUp.m` - calculates the average heat content in the upper 40 m at modeled freeze up for profiles inside and outside the meltwater - these values are used in the discussion section of the paper.

`calculateCompassAngle.m` - Takes in zonal and meridional velocity vector components and calculates the compass angle of that velocity (i.e. north = 0°, east = 90°). Used in `plot_frontTracking.m`

`calculateHeatContent.m` - Depth-integrates the heat content relative to freezing of a T-S profile. Uses salinity data to calculate freezing point. If surface data is missing, fills with the shallowest measurement. Then if additional data is missing fills with linear interpolation. 

`defineSODAconstants.m` - Various useful values defined in one place: standard colors used throughout, locations of the moorings within the study area
