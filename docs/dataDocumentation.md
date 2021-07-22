## Observational data and PWP model results 

This document describes the observational data (from Seaglider, underway CTD, Wave Glider, and USCGC *Healy's* underway system) and PWP model results used in “Meltwater advection hastens autumn freeze up” by Crews et al. These data are available on the UW ResearchWorks archive at http://hdl.handle.net/1773/47135

Conservative temperature and absolute salinity were calculated using the Gibbs SeaWater (GSW) Oceanographic Toolbox (v3.0.11) functions `gsw_CT_from_t.m` and `gsw_SA_Sstar_from_SP.m` available at http://www.TEOS-10.org. 

Data are included in the following files (last updated 3 July, 2021 by Laura Crews lcrews@uw.edu).

* healyUnderway.csv
* waveglider.nc
* profiles.nc
* pwpResults_profile*.nc (multiple files with profile number *)

Descriptions of each of the data files are given below.

### healyUnderway.nc
This file contains data collected by USCGC *Healy’s* seawater intake system while underway. It includes raw ocean temperature and salinity data as well as ocean temperature and salinity after the corrections described in Crews et al. (2021) section 2.2.3 have been applied. Conservative temperature and absolute salinity were calculated first, then the corrections were calculated and applied. Note that there are a few clearly erroneous high temperature and low salinity spikes in the underway data that were excluded from our analysis and should be removed by future users of this data. 

Variables:
* lat: Latitude [degrees N]
* lon: Longitude [degrees E]
* time: Time of measurement in Matlab date format [days from Jan 0, 0000]
* temps_uncorrected: Raw temperature data from the underway system [degrees C]
* salts_uncorrected: Raw practical salinity data from the underway system [psu]
* CTs_uncorrected: Conservative temperature calculated from the uncorrected absolute salinity and uncorrected raw temperature [degrees C]
* SAs_uncorrected: Absolute salinity calculated from the uncorrected practical salinity [g/kg]
* temps: Corrected temperature data (raw in situ temperature with the median difference between the uCTD temperature and co-located underway temperature added) [degrees C]
* salts: Corrected practical salinity data (raw practical salinity with the correction derived from a linear fit of the differences between the uCTD practical salinity and co-located underway practical salinity added) [psu]
* CTs: Corrected conservative temperature data (CTs_uncorrected with the median difference between the uCTD conservative temperature and co-located uncorrected underway conservative temperature added) [degrees C]
* SAs: Corrected absolute salinity data (SAs_uncorrected with the correction derived from a linear fit of the differences between the uCTD absolute salinity and co-located underway absolute salinity added) [psu]

### waveglider.nc
This file contains ocean temperature and salinity data collected by four Wave Glider autonomous vehicles. There are three potential measurement depths given in the “z” variable; not all vehicles collected data at all depths. Conservative temperature and absolute salinity are calculated at the depths for which conductivity data are available. 

* vehicle: Wave Glider ID number
* lats: Latitude [degrees N]
* lons: Longitude [degrees E]
* times: Time of measurement in Matlab date format [days from Jan 0, 0000]
* z: Depths of measurements (3 possible depths per vehicle) below ocean surface (positive downwards) [m]
* temps: In situ water temperature [degrees C]
* salts: In situ practical salinity [psu]
* CT: Conservative temperature [degrees C]
* SA: Absolute salinity [g/kg]
* sigthes: Potential density referenced to the surface [kg/m3]

### profiles.nc
This file contains the Seaglider profiles and underway CTD profiles sorted by measurement time. The “vehicle” field indicates which Seaglider took each profile, or if that profile was taken by the uCTD. Profiles lacking data shallower than 15 m are not included. Missing surface data have not been filled, but were filled as a constant value for the analyses in Crews et al. (2021). The “qualFlag” field indicates if a profile was excluded from the Crews et al. (2021) analyses due to density inversions in the upper 40 m. Note that profiles may have poor data quality at deeper depths. No filter was applied to these depths as the present study was concerned with the surface ocean, but future users should inspect profiles for problems at depth. 

* vehicle: Seaglider ID# or, if equal to 1, indicates uCTD
* qualFlag: This flag value = 1 if a profile was considered "good" or this flag = 0 if this profile was excluded from analysis. Profiles were filtered as “good” if they lacked density inversions larger than 0.05 kg/m3 in the upper 40 m (the depth range important for the analysis done in this paper). 
* lats: Latitude [degrees N]
* lons: Longitude [degrees E]
* times: Time of measurement in Matlab date format [days from Jan 0, 0000]
* z: Depth below ocean surface (positive downward) [m]
* temps: Temperature [degrees C]
* salts: Practical salinity [psu]
* CT: Conservative temperature [degrees C]
* SA: Absolute salinity [g/kg]
* sigthes: Potential density referenced to the surface [kg/m3]
* pwpFreezeTime: Time the 3rd depth bin (centered at 2.5 m) of the PWP model initiated with this profile fell below the freezing temperature, as a proxy for freeze up [days from Jan 0, 0000] 
* pwpz: Center depth of PWP model depth bin (positive downward) [m]
* pwpCT: PWP modeled conservative temperature at pwpFreezeTime [degrees C]
* pwpSA: PWP modeled absolute salinity at pwpFreezeTime [g/kg]
* pwpSigthes: PWP modeled potential density at pwpFreezeTime [kg/m3]

### pwpResults_profile*.nc
These files contain the temperature and salinity output for the PWP model initiated with the observed profile number indicated in the file name. The observed profile number corresponds to the “pwpProfNum” field in the profiles.nc file. 

Each file contains the **global attributes**:
* profile Number: the profile number also indicated in the file name which maps to the data in the profiles.nc file of the same “pwpProfNum”
* lat, lon, and time of observed profile: Latitude, longitude, and time of the Seaglider or uCTD profile used to initiate this PWP model run, copied from the profiles.nc file for this profile number
* heatBudgetTransect: Indicates which transect this profile was on for the heat budget (i.e., which row of Figures 8 and 9 of Crews et al. (submitted 2021) these data were used in). 

Each file also contains the  **variables**:
* time: Time of PWP model timestep [Days from Jan 0, 0000]
* z: Center depth of PWP model depth bin (positive downward) [m]
* heatFlux: Net surface heat flux from ERA5 (sum or net shortwave, net longwave, sensible, and latent components) at each model time step [W/m2]
* CT: Conservative temperature of model output [degrees C]
* SA: Absolute salinity of model output [g/kg]
* sigthe: Potential density of model output referenced to the surface [kg/m3]
