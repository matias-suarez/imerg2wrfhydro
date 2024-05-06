# imerg2wrfhydro

Scripts for regridding IMERG precipitation to the WRF-Hydro domain.
This script uses inverse distance interpolation. A precipitation value is calculated for each WRF-Hydro grid point using the 4 nearest imerg grid points.
Files required:
1) IMERG files in netcdf format
2) geogrid file (this file contains the WRF-Hydro domain information)

Note: This repository contains an .ipynb to download the IMERG precipitation from the GES DISC website (https://disc.gsfc.nasa.gov/).
