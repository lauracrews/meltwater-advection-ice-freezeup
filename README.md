# Meltwater advection hastens autumn freeze up
Code to reproduce analysis and figures of "Meltwater Advection Hastens Freeze Up" by Crews et al. (2021)

# Dependencies

### MATLAB and toolboxes
* MATLAB (R2018b used to develop this code) <br />
* MATLAB image processing toolbox function imgaussfilt.m <br />
* Gibbs SeaWater (GSW) Oceanographic Toolbox Version 3.0.11, available at http://www.TEOS-10.org <br />
* m_map toolbox, available at https://www.eoas.ubc.ca/~rich/map.html#download <br />

### Colormaps used in figures
* cmocean, available on the MATLAB file exchange at https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps <br />
* ColorBrewer, available on the MATLAB file exchange at https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps <br />

# External data

Download the **sea surface temperature** data from the Physical Oceanography Distributed Active Archive Center at NASA/JPL. You will need to make a free account. <br />

For 8-day files go to: 
https://podaac-tools.jpl.nasa.gov/drive/files/allData/modis/L3/terra/11um/v2019.0/4km/8day/2018 <br />

Download the files: <br />
T20182572018264.L3m_8D_NSST_sst_4km.nc <br />
T20182732018280.L3m_8D_NSST_sst_4km.nc <br />
T20182652018272.L3m_8D_NSST_sst_4km.nc <br />
T20182812018288.L3m_8D_NSST_sst_4km.nc <br />

For daily files go to:
https://podaac-tools.jpl.nasa.gov/drive/files/allData/modis/L3/terra/11um/v2019.0/4km/daily/2018  <br />

Download the files for 19 September 2018 to 16 October 2018:  <br />
TERRA_MODIS.20180919.L3m.DAY.NSST.sst.4km.nc  <br />
…  <br />
TERRA_MODIS.20181016.L3m.DAY.NSST.sst.4km.nc  <br />

Download the **sea surface salinity** data from https://www.seanoe.org/data/00607/71909/ <br />
Scroll down and click the Weekly SMOS Arctic SSS v1.1 link to download a large .zip file, from which you should extract the files:

SMOS-arctic-LOCEAN-SSS-2018-09-15-v1.1AT-7days.nc	<br />
SMOS-arctic-LOCEAN-SSS-2018-09-21-v1.1AT-7days.nc	<br />
SMOS-arctic-LOCEAN-SSS-2018-10-07-v1.1AT-7days.nc <br />
SMOS-arctic-LOCEAN-SSS-2018-09-30-v1.1AT-7days.nc <br />

Download the CryoSAT2 **Dynamic ocean topography** data from http://rads.tudelft.nl/rads/data/authentication.cgi <br />
Navigate through three pages of options and make the following selections, organized by page:
1. You will need to enter your email address to receive the data. Use the default settings under "Options" <br />
1. Select the "sea level anomaly" variable <br />
1. The date range of interest is covered by Cycles 109 & 110. Select all passes. Update the geographic range <br />


