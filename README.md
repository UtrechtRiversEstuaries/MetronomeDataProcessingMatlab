## MetronomeDataProcessingMatlab
Files for processing raw data from the tidal tilting flume 'The Metronome' in Matlab

# DEM processing
* Laser2NetCDF.m
This file processes the raw laser scanner .tif images to DEMs in netCDF format. All DEMs of one experiment are saved in a single file including metadata. The file is adjusted to the standard folder structure.

* laserMetronome_ESL.m
This file is a function used by Laser2NetCDF.m. Here the actual processing takes place.

* NetCDF2Plots
This file creates basic plots from Metronome DEMs in NetCDF format.

