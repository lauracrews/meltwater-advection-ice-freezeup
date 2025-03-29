## Overview
This directory contains code used to process the observational data from the shipboard sensors, underway CTD, Seagliders, and Wave Gliders. 
It produces the data products archived at [UW ResearchWorks](http://hdl.handle.net/1773/47135) and described in the [observational data documentation](https://github.com/lauracrews/meltwaterAdvection/blob/main/docs/dataDocumentation.md).

This code is provided for reference and completeness. **It is not necessary to rerun this code - instead, the processed data should be downloaded from the archive at [UW ResearchWorks](http://hdl.handle.net/1773/47135).** 

The master script for data processsing is `process_observational_data.m`. It calls subroutines necessary for data processsing and saves the datasets later archived on ResearchWorks. Tasks performed:
* Corrects the temperature and salinity data from Healy's underway system, writes the file `healyUnderway_test.csv`
* Loads Wave Glider data, calculates conservative temperature and absolute salinity, writes the file `waveglider.nc`
* Creates a combined data structure of all uCTD an Seaglider profiles in the specified region of interest, calculates the heat content of these, flags low quality data, writes the file `profiles.nc`
* Initiates the PWP model with observed temperature and salinity profiles. Runs the PWP model with ERA5 atmospheric forcing near the observation used to initialize the model. Saves individual .nc files for each PWP model run used in the head budget 

## Wave Gliders
`loadWaveglider.m` - Iterates through all Wave Glider files and combines the desired datasets (ocean temperature and salinity) into one data structure. Calculates conservative temperature and absolute salinity. Can be updated to also return the Wave Glider meteorological data. 

## Seaglider and uCTD profiles
`profilesInRegion_timeseries.m` - Makes a data structure of temperature and salinity profiles collected by all platforms (all Seagliders and uCTD) in a specified geographic region and date range. Calls `calculateHeatContent.m` to depth-integrate the heat in the temperature profiles. This is run once to identify all of the profiles used in the study, after which the data structure is saved and the relevant PWP model results are added later

## PWP modeling from observed profiles
`initiate_run_pwp.m` - Prepares to run the PWP model. Averages ERA5 meteorological variables near each set of observed profiles used to initialize the model. Runs model. Saves PWP model output at each timestep for each observational profile, as well as adding key PWP results (i.e. profiles at freezing and freeze up time) to the `profiles.mat` data structure. The resulting intermediate data product is `profiles_pwpresults_withFreshwater.mat', available in the ~/data/process_raw_observational_data/ directory of this repo

`makeERA5timeseries.m` - Averages ERA5 data within the specified geographic coordinates. Input the names of the desired parameters as strings

`pwp_timeseries.m` - Initializes and runs the PWP model code. Intakes data structure of atmospheric data with fields for surface stress, heat fluxes, precip/evap. Initializes PWP with observed temperature and salinity profiles. Sets various required parameters (diffusivities, time stepping, LW and SW extinction coefficients). After running model, calculates heat content of modeled profiles. 

## Underway meteorological and near-surface temperature and salinity from USCGC Healy
`load_correct_underwayData.m` - Loads temperature and salinity data collected by the ship's seawater intake system while underway. Compares Healy underway temperature and salinity data to concurrent uCTD data and applies a linear fit correction to the underway salinity data and a constant correction to the underway temperature data (see Alory et al., 2015, [doi:10.1016/j.dsr.2015.08.005](https://www.sciencedirect.com/science/article/pii/S0967063715001417)). These data are saved as `healyUnderway.csv`

`load_HEALY1802_MET.m` - Luc Rainville’s script which implements San Nguyen’s script `SN_readShipMET.m` to extract desired data fields at one minute intervals from the raw Healy underway data which was provided to Laura Crews by Luc Rainville. The resulting .mat file is `met.mat` Extracted data fields are summarized in the table below and descriptions of all data fields are included in `SN_readShipMET.m.` In addition to meteorological data, which were not used in the analysis, this script **extracts the underway temperature and salinity data** which were used in the analysis

Note that for some meteorological parameters Healy has [multiple sensors](https://github.com/lauracrews/meltwaterAdvection/blob/main/code/process_raw_observational_data/Shipboard%20Science%20Data%20Collection%20Map_HCO%20and%20Transducer%20View.pdf) and I have not been able to determine which sensor was used for each reported measurement. This would be important for determining if wind measurements were sheltered by the ship’s structure and for using the correct sensor heights when calculating turbulent heat fluxes using COARE. I have not pursued this further because I use data from ERA5 in the subsequent analyses as I need meteorological data even when Healy was not located in the study area. However the problem of the unknown sensor locations and heights should be considered further if Healy’s meteorological data is used for analysis in the future. 

| Variable Description from `SN_readShipMET.m`  | Field Name                                          | Notes |
|-------------|--------------------------------------------------------|-------|
| LA          | NMEA Latitude                                          |       |
| LO          | NMEA Longitude                                         |       |
| TT          | TSG Temperature (°C)                                   | Temperature of water from seawater intake, located 2.7 m below Healy’s design waterline |
| SA          | Salinity (psu)                                         | Salinity of water from seawater intake |
| AT          | Air Temperature (°C)                                   |       |
| BP          | Barometric Pressure (mb)                               |       |
| TW          | True Wind Speed (m/s)                                  |       |
| TI          | True Wind Direction                                    | Direction wind is coming from |
| LW          | Long wave radiation (W/m²)                             |       |
| SW          | Short wave radiation (W/m²)                            |       |
| RH          | Relative Humidity (%rh)                                |       |


