## MetronomeDataProcessingMatlab
Files for processing raw data from the tidal tilting flume 'The Metronome' in Matlab

# DEM processing
* Laser2NetCDF.m
This file processes the raw laser scanner .tif images to DEMs in netCDF format. All DEMs of one experiment are saved in a single file including metadata. The file is adjusted to the standard folder structure.

* laserMetronome_ESL.m
This file is a function used by Laser2NetCDF.m. Here the actual processing takes place.

* NetCDF2Plots.m
This file creates basic plots from Metronome DEMs in NetCDF format.

# Laser gantry camera processing
* Stitch_SLRcameraLaserGantry_ESL_v0.m
This file combines the photos of the SLR camera on the laser gantry into one photo. It is still adjusted for experiments with half the flume and the camera positions and lighting differences are not yet perfectly calibrated. It also is not adjusted for the standard folder structure.

# Overhead cameras processing
* Stitch_overheadcams_ESL_v0.m
This file combines the photos of the overhead cameras into one photo and also processes them to videos. It is still adjusted for experiments with half the flume and the camera positions and lighting differences are not yet perfectly calibrated. It also is not adjusted for the standard folder structure.

# Water level measurements processing
No files for this are currently included as the files that were used in previous experiments are quite messy and not well adjusted for general processing. A better version will be coded in Python.
