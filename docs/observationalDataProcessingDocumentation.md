This file describes code used to process the observational data from the shipboard instrumentation, underway CTD, Seagliders, and Wavegliders. This code is provided for reference and completeness.

`process_observational_data.m` - Calls subroutines necessary for data processsing

**It is not necessary to rerun this code - instead, the processed data should be downloaded from the archive at [UW ResearchWorks](http://hdl.handle.net/1773/47135).** 

## Wave Gliders
`loadWaveglider.m` - Iterates through all Wave Glider files and combines the desired datasets (ocean temperature and salinity) into one data structure. Calculates conserevative temperature and absolute salinity. Can be updated to also return the Wave Glider meteorological data. 

## Seaglider and uCTD profiles
`profilesInRegion_timeseries.m` - Makes a data structure of temperature and salinity profiles collected by all platforms (all Seagliders and uCTD) in a specified geographic region and date range. Calls calculateHeatContent to depth-integrate the heat in the temperature profiles. This is run once to identify all of the profiles used in the study, after which the data structure is saved and the relevant PWP model results are added later

`loadProfiles.m` -  Loads the data structure of all Seaglider and uCTD profiles in the study. Fills missing surface data and eliminates profiles with density inversions

`identifyProfileNumbers.m` - Function to identify Seaglider and uCTD profiles in specified geographic bounds and time period

## Underway meteorological and near-surface temperature and salinity from USCGC Healy
`load_correct_underwayData.m` - Loads temperature and salinity data collected by the ship's seawater intake system while underway. Compares Healy underway temperature and salinity data to concurrent uCTD data and applies a linear fit correction to the underway salinity data and a constant correction to the underway temperature data (see Alory et al., 2015, [doi:10.1016/j.dsr.2015.08.005](https://www.sciencedirect.com/science/article/pii/S0967063715001417)). These data are saved as 'healyUnderway.csv'

`load_HEALY1802_MET.m` - Luc Rainville’s script which implements San Nguyen’s script `SN_readShipMET.m` to extract desired data fields at one minute intervals from the raw Healy underway data which was provided to Laura Crews by Luc Rainville. The resulting .mat file is `met.mat` Extracted data fields are summarized in the table below and descriptions of all data fields are included in `SN_readShipMET.m.`

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


