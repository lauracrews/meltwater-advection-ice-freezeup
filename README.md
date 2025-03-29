# Meltwater advection hastens autumn freeze up
This repo contains code to reproduce the analysis and figures in the paper "Direct observations of the role of lateral advection of sea ice meltwater in the onset of autumn freeze up" by Crews et al. (2022)

Crews, L., Lee, C. M., Rainville, L., & Thomson, J. (2022). Direct observations of the role of lateral advection of sea ice meltwater in the onset of autumn freeze up. Journal of Geophysical Research: Oceans, 127, e2021JC017775. https://doi.org/10.1029/2021JC017775

### Key Points
* High spatial and temporal resolution observations by ship-based and autonomous instruments in the ice-free Beaufort Sea are presented
* Modified meltwater advected about 100 km over several weeks, cooling and shoaling the mixed layer and hastening freeze up by several days
* Meltwater advection caused nearly as much mixed layer heat loss as was caused by seasonally integrated heat loss to the atmosphere

### Plain Language Summary
In large parts of the Arctic Ocean, sea ice melts completely in summer and reforms in autumn in a process known as “freeze up.” The timing of freeze up may affect energy exchanges between the ocean and atmosphere, such as occur during storms, as well as impact the sea ice ecosystem. Warming ocean temperatures mean freeze up is occurring later, lengthening the season in which vessels that are not ice-rated can engage in operations like shipping, resource extraction, and supplying local communities with fuel and cargo. Accurate freeze up forecasting helps ensure these operations are conducted safely. This study uses ocean temperature and salinity data from autonomous vehicles as well as data from satellites to demonstrate how areas that were affected by recent sea ice melt can freeze up early. We observed a cold and fresh footprint of meltwater near sea ice in the Beaufort Sea that flowed away from its original location and altered the upper ocean in neighboring areas. A simple model that simulates the atmosphere's influence on the ocean demonstrates that this meltwater caused the ocean to freeze several days earlier than it otherwise would have by cooling the upper ocean and limiting the depth of ocean mixing.
![ice_formation_melt](https://github.com/lauracrews/meltwaterAdvection/blob/main/docs/fig1/ice_formation_melt.png)

# Dependencies

### MATLAB and toolboxes
* MATLAB (R2018b used to develop this code) <br />
* MATLAB image processing toolbox (the function imgaussfilt.m is used) <br />
* Gibbs SeaWater (GSW) Oceanographic Toolbox (version 3.0.11 used here), available at http://www.TEOS-10.org <br />
* m_map Toolbox (version 1.4m used here), available at https://www.eoas.ubc.ca/~rich/map.html#download <br />

### Colormaps used in figures
* cmocean, available on the MATLAB file exchange at https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps <br />
* ColorBrewer, available on the MATLAB file exchange at https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps <br />

# Getting started
Reproducing the analysis and figures requires the following steps:
1. Download this repository to the first directory on the Matlab search path (this can be found by running `userpath`). <br /> **Note:** If you would like to store this repo elsewhere, edit the variable `rootPath` in the script `run_meltwaterAdvection.m` to reflect the location of this repository’s contents
2. Download the needed external data to the appropriate directories (instructions below)
3. Within Matlab, add the ~/meltwaterAdvection/ directory and subdirectories to the path
4. Enter `run_meltwaterAdvection` in the Matlab command line to reproduce all figures 

## Download data
From the University of Washington ResearchWorks archive at http://hdl.handle.net/1773/47135, download 
* Seaglider, underway CTD, Wave Glider, and USCGC *Healy* underway data to the directory ~/meltwaterAdvection/data/ 
* PWP model results used to make the heat budget should be unzipped into the directory ~/meltwaterAdvection/data/PWPresults/ 

Descriptions of these data are available in the [observational data documentation](https://github.com/lauracrews/meltwaterAdvection/blob/main/docs/dataDocumentation.md). The processing routines for the raw observational data and PWP model are described in the [observational data processing documentation](https://github.com/lauracrews/meltwater-advection-ice-freezeup/blob/main/code/process_raw_observational_data/README.md).

**It is not necessary to rerun the observational data processing and PWP model code - instead, the processed data should be downloaded from the archive at [UW ResearchWorks](http://hdl.handle.net/1773/47135).** 

##
**Sea ice concentration** from AMSR2 will be downloaded from the University of Bremen [archive](https://seaice.uni-bremen.de/data/amsr2/asi_daygrid_swath/n3125/2018/)  by `batchAMSR2` within `run_meltwaterAdvection`  <br />
##
Download the MODIS-Terra **sea surface temperature** data from the Physical Oceanography Distributed Active Archive Center at NASA/JPL. You will need to make a free account. <br />

For 8-day files go to: 
https://podaac-tools.jpl.nasa.gov/drive/files/allData/modis/L3/terra/11um/v2019.0/4km/8day/2018 <br />

Download the following files to the directory ~/meltwaterAdvection/data/ModisTerra/ <br />
T20182572018264.L3m_8D_NSST_sst_4km.nc <br />
T20182732018280.L3m_8D_NSST_sst_4km.nc <br />
T20182652018272.L3m_8D_NSST_sst_4km.nc <br />
T20182812018288.L3m_8D_NSST_sst_4km.nc <br />

For daily files go to:
https://podaac-tools.jpl.nasa.gov/drive/files/allData/modis/L3/terra/11um/v2019.0/4km/daily/2018  <br />

Download the files for 19 September 2018 to 16 October 2018 to the directory ~/meltwaterAdvection/data/ModisTerra/daily/  <br />
TERRA_MODIS.20180919.L3m.DAY.NSST.sst.4km.nc  <br />
…  <br />
TERRA_MODIS.20181016.L3m.DAY.NSST.sst.4km.nc  <br />
##
Download the SMOS **sea surface salinity** data from https://www.seanoe.org/data/00607/71909/ <br />
Scroll down and click the Weekly SMOS Arctic SSS v1.1 link to download a large .zip file, from which you should extract the following files to the directory ~/meltwaterAdvection/data/SMOS_SSS/

SMOS-arctic-LOCEAN-SSS-2018-09-15-v1.1AT-7days.nc	<br />
SMOS-arctic-LOCEAN-SSS-2018-09-21-v1.1AT-7days.nc	<br />
SMOS-arctic-LOCEAN-SSS-2018-09-30-v1.1AT-7days.nc <br />
SMOS-arctic-LOCEAN-SSS-2018-10-07-v1.1AT-7days.nc <br />
##
Download the CryoSat-2 **dynamic ocean topography (DOT)** data from the Radar Altimeter Database System at http://rads.tudelft.nl/rads/data/authentication.cgi <br />
You will need to enter your email address to receive the data. Navigate through three pages of options and make the following selections, organized by page:
1. Use the default settings under "Options" <br />
1. Select the "sea level anomaly" variable <br />
1. The date range of interest is covered by Cycles 109 & 110. Select all passes. Update the geographic range to 70°N to 90°N, -180°E to -130°E <br />

You will then receive an email to download several hundred zipped .asc files. 

**Alternatively, the DOT data used in the study are available on my github as [rawDOT.zip](https://github.com/lauracrews/meltwaterAdvection/raw/main/rawDOT.zip)**

Unzip the files obtained by either method into the directory ~/meltwaterAdvection/data/DOT/raw/

These data will be read into a file called `allDOT.mat` in the course of the analysis. A copy of my version of this file is in the data directory of this repo.  

## Run analysis and create figures

The script `run_meltwaterAdvection.m` will create all the figures in the paper and save them as .png and Matlab .fig files into their own subdirectories within ~/meltwaterAdvection/figures/ 

**Note:** Copies of the output figures are available in numbered subdirectories in the [docs](https://github.com/lauracrews/meltwaterAdvection/tree/main/docs) directory of this repository if you would like to preview the results. 

Switches determining if figures should be saved are provided within each plotting script (`plot_*.m`). If you do not want to save the figures, edit the plotting scripts directly and set `saveFigs = false`

Descriptions of each of the subroutines called by `run_meltwaterAdvection.m` are in the [code documentation](https://github.com/lauracrews/meltwaterAdvection/blob/main/docs/codeDocumentation.md)
