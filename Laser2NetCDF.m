% This script processes the raw .tif-files from the Mtronome laser scanner
% into a single netCDF-file per experiment/pilot and automatically appends
% the processed DEM if such a file with DEMs of previous timesteps already
% exists. 
% You have to specify your personal working directory, the number of the
% current experiment, your desired resolution and range of interest,
% whether you want to save the DEM as .mat-file as well and all the
% attributes containing the metadata about the experiment and the DEM. You
% also have to make sure that the folders in the laser_scanner directory
% withthe raw data have names only containing the cycle number but nothing
% else. In each folder there should be exactly 20 files, otherwise an error
% will occur. So make sure to get rid of the extra files the laser scanner
% always saves (usually the last one).
% Jan-Eike Rossius, Utrecht University in October 2022

close all
clear

%% Input - to be filled in individually
% your personal working directory in which you use the standard folder structure:
pwd = 'C:\Users\7062052\OneDrive - Universiteit Utrecht\Open Science Metronome\MetronomeWorkingDirectory\';
% the current experiment (three digit number, in case of pilots use e.g. '052\Pilot1' as this string is used for directory paths)
exp = '052\Pilot4';
% the desired grid resolution of the processed DEM in metres:
gridresol = 3/1000; 
% range of interes within the metric flume coordinates in the form [xmin xmax ymin ymax]:
roi = [0 20 0 3]; 
% the percentiles for the low and high data of the DEM:
percen = [10 90];
% Option to save the individual DEMs also as .mat-file
matsave = true;

% --------------- %  Attributes  % --------------- %

% Global attributes
title='Data of fourth pilot for flood-dominant equilibrium estuaries';
long_title='Data from different sources of the fourth pilot for flood-dominant equilibrium estuaries in the Metronome';
comments='File in development';
institution='Utrecht University';
source='Metronome instruments';
references='none';
history='laser data processed to DEMs and written to NetCDF';
summary='to be filled in';
keywords='flume experiments, physical modelling, tidal channels, estuaries, tidal assymmetry, steady state';

conventions='mostly CF-1.8';  
conventionshelp='http://cfconvention.org/Data/cf-conventions/cf-conventions-1.9/cf-conventions.html';

author='Jan-Eike Rossius';
correspondence='Prof. Maarten G. Kleinhans, m.g.kleinhans@uu.nl';

% DEM group attributes
dems_title='DEMs of fourth pilot for flood-dominant equilibrium estuaries';
dems_long_title='Digital elevation models of the fourth pilot for flood-dominant equilibrium estuaries in the Metronome';
dems_history='laser data processed to DEMs and written to NetCDF';
% dems_history='laser data processed to DEMs and written to NetCDF (the data of 3001 cycles is the scan that was amde with water in the system)';
dems_description='This data contains the elevation data of the sediment bed in the tidal flume The Metronome. The raw data was measured with a laser scanner which was then processed to give calibrated elevation data above the flume floor. Since several points (usually 9) are processed into one grid cell, the median and 10th and 90th percentile are given.';

% Metronome attributes
% General
expnr='052_P4';
labslaves='Eise Nota and Jan-Eike Rossius';
metronome_comments='long channel, no barriers, flood asym, no feed';
manmade_structures='no';
nr_of_cycles=5000;
start_of_experiment=datestr(2022/10/13);
end_of_experiment=datestr(2022/10/24);
% Boundary conditions - tilting and weir
amplitude_tilt_1='75 mm';
period_tilt_1='40 s';
phase_tilt_1='0 degrees';
offset_tilt_1='0 mm';
amplitude_tilt_2_overtide='15 mm';
period_tilt_2_overtide='20 s';
phase_tilt_2_overtide='90 degrees';
amplitude_weir='12 mm';
period_weir='40 s';
offset_weir='65 mm';
phase_weir_in_relation_to_tilting='180 degrees';
amplitude_weir_overtide='2.4 mm';
period_weir_overtide='20 s';
phase_weir_overtide='90 degrees';
flood_ebb_dominance='flood';
% Boundary conditions - river
river_water_discharge='0 l/h';
river_water_timing='none';
river_sand_feed_rate='0 l/h';
river_sand_timing='none';
river_walnut_feed_rate='0 l/h';
river_walnut_timing='none';
% Boundary conditions - waves
waves_on_off='on';
waves_frequency='2 Hz';
waves_timing='flood only';
% Boundary conditions - sediment at inlet
sediment_at_inlet_yes_no='no';
inlet_sand_feed_rate='0 l/h';
inlet_sand_timing='none';
inlet_walnut_feed_rate='0 l/h';
inlet_walnut_timing='none';
% Boundary conditions - sea-level rise
sea_level_rise_yes_no='no';
sea_level_rise_start_cycle='none';
sea_level_rise_end_cycle='none';
sea_level_rise_cycle_increment='none';
sea_level_rise_weir_change_increment='none';
sea_level_rise_total_nr_of_events=0;
% Initial conditions
initial_sediment_bed_level='100 mm';
initial_coastline_x_position='18 m';
initial_channel_depth='30 mm';
initial_channel_width='540 mm';
initial_channel_width_variation='none';
initial_tidal_basin_depth='none';
initial_tidal_basin_width='none';
initial_tidal_basin_length='none';
initial_barriers_present='no';
initial_barrier_width='none';
initial_barrier_height='none';
initial_inlet_width='none';
initial_barriers_fixated='no';
initial_landward_barrier_present='no';
initial_landward_barrier_x_position='none';
initial_landward_barrier_fixated='no';
% Vegetation
vegetation_yes_no='no';
vegetation_start_sowing='none';
vegetation_interval_sowing='none';
vegetation_nr_of_sowing_events=0;
vegetation_river_Lotus='0 g';
vegetation_river_Alfalfa='0 g';
vegetation_river_Veronica='0 g';
vegetation_inlet_Lotus='0 g';
vegetation_inlet_Alfalfa='0 g';
vegetation_inlet_Veronica='0 g';

%% Evaluating the directory
folders=dir([pwd '01metronome_experiments\Exp' exp '\raw_data\laser_scanner\']); % find all subfolders in the laser data directory
if isempty(folders) % throw error if there are none
    error("Directory path incorrect or no data available")
elseif sum(vertcat(folders.isdir))<length(folders) % Throw error if there are files not in subfolders
    error("There are files not in subfolders per cycle. Make sure they are in the correct subfolder")
end

% Find the elements . and .. in the folders
nofolderindices=0;
for i=1:length(folders)
    idx=find(folders(i).name=='.');
    if ~isempty(idx)
        if nofolderindices==0
            nofolderindices=i;
        else
            nofolderindices=[nofolderindices i];
        end
    end
end

% Get rid of . and ..
if nofolderindices(1)>1
    for i=1:nofolderindices(1)-1
        foldersnew(i)=folders(i);
    end
end
for i=nofolderindices(2)+1:length(folders)
    foldersnew(i-2)=folders(i);
end

%% Definintion of important variables
expalt=strrep(exp,'\','-'); % repacing \ in exp string to have a string that is nicer for naming files etc.

% Extracting numerical timesteps from folder names
timesteps=NaN(length(foldersnew),1);
for i=1:length(foldersnew)
    timesteps(i)=str2double(foldersnew(i).name);
end

% defining x and y coordinates from range of interest
x=[roi(1):gridresol:roi(2)];
nx=length(x)-1;
y=[roi(3):gridresol:roi(4)];
ny=length(y)-1;

% Value for no data
fillValue = -9999;

%% Checking existing file or creating new one
try % Check whether file in question exists and if so, open it
    fileid=netcdf.open([pwd '01metronome_experiments\Exp' exp '\processed_data\DEMs\laser_scanner\Exp' expalt 'DEMs.nc'],'WRITE');
    new=false; % variable if new file was created for later
catch % create new file if it does not exist yet
    fileid=netcdf.create([pwd '01metronome_experiments\Exp' exp '\processed_data\DEMs\laser_scanner\Exp' expalt 'DEMs.nc'],'NETCDF4');
    new=true; % variable if new file was created for later
end

%% Basic file information and dimensions
if new==true % in new file, define information and dimensions and fill them with data
    % Define group
    dems=netcdf.defGrp(fileid,'DEMs');
    
    % Define useful constants:
    NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');

    % Define dimensions:
    dimidT = netcdf.defDim(dems,'time',netcdf.getConstant('NC_UNLIMITED'));
    dimidY = netcdf.defDim(dems,'y',ny);
    dimidX = netcdf.defDim(dems,'x',nx);
    dimidN = netcdf.defDim(dems,'nv',2);
    
    % Define axis:
    xdim = netcdf.defVar(dems,'x','double',[dimidX]);         
    netcdf.putAtt(dems,xdim,'long_name','distance along the flume');
    % netcdf.putAtt(dems,xaxis,'standard_name','longitude');
    netcdf.putAtt(dems,xdim,'units','m');
    netcdf.putAtt(dems,xdim,'valid_range',[0 20]);
    netcdf.putAtt(dems,xdim,'axis','X');
    netcdf.putAtt(dems,xdim,'bounds','x_bounds');
    netcdf.defVarFill(dems,xdim,false,fillValue);
    netcdf.endDef(fileid);
    netcdf.putVar(dems,xdim,x(1:end-1)+gridresol/2);
    
    ydim = netcdf.defVar(dems,'y','double',[dimidY]);         
    netcdf.putAtt(dems,ydim,'long_name','distance across the flume');
    % netcdf.putAtt(dems,yaxis,'standard_name','latitude');
    netcdf.putAtt(dems,ydim,'units','m');
    netcdf.putAtt(dems,ydim,'valid_range',[0 3]);
    netcdf.putAtt(dems,ydim,'axis','Y');
    netcdf.putAtt(dems,ydim,'bounds','y_bounds');
    netcdf.defVarFill(dems,ydim,false,fillValue);
    netcdf.endDef(fileid);
    netcdf.putVar(dems,ydim,y(1:end-1)+gridresol/2);
    
    tdim = netcdf.defVar(dems,'time','double',[dimidT]);         
    netcdf.putAtt(dems,tdim,'long_name','time since the start of the experiment');
    netcdf.putAtt(dems,tdim,'units','tidal_cycles');
    netcdf.putAtt(dems,tdim,'valid_range',[0 20000]);
    netcdf.putAtt(dems,tdim,'axis','T');
    netcdf.defVarFill(dems,tdim,false,fillValue);
    netcdf.endDef(fileid);
%     netcdf.putVar(dems,taxis,timesteps);
    
    %Defining boundary variables
    xbnds = netcdf.defVar(dems,'x_bounds','double',[dimidN, dimidX]);         
    netcdf.defVarFill(dems,xbnds,false,fillValue);
    netcdf.putVar(dems,xbnds,[x(1:end-1)' x(2:end)']);
    
    ybnds = netcdf.defVar(dems,'y_bounds','double',[dimidN, dimidY]);         
    netcdf.defVarFill(dems,ybnds,false,fillValue);
    netcdf.putVar(dems,ybnds,[y(1:end-1)' y(2:end)']);

    % -----------------------------------------------%

else % if file already exists, extract information, dimensions and data

    % Get information about the contents of the file.
    [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(fileid);
    
    % Identify location where DEMs are stored
    if numdims==0 % If parent group of file doesn't contain any dimensions, search in groups
        childGrps = netcdf.inqGrps(fileid);
        for i=1:length(childGrps) % loop through groups
            groupName = netcdf.inqGrpName(childGrps(i));
            if strcmp(groupName,'DEMs') % find group called 'DEMs'
                dems=childGrps(i);
                [numDEMdims, numDEMvars, numDEMglobalatts, unlimDEMdimID] = netcdf.inq(dems); % get info about the DEM group
            end
        end
    else % if no DEM group exists, the parent group acts as DEM group
        dems=fileid;
        [numDEMdims, numDEMvars, numDEMglobalatts, unlimDEMdimID] = netcdf.inq(dems);
    end

    for j=0:numDEMdims-1 % loop through dimensions to identify which dimension is which
        try
            axisAtt = netcdf.getAtt(dems,j,'axis'); % get the axis attribute
        catch exception
            axisAtt = 'none'; % put 'none' if axis attribute doesn't exist
        end
        if strcmp(axisAtt,'X') % if axis attribute is 'X', the current dimension is the x-dimension
            xdim=j; % define the dimension id
            xdimdata=netcdf.getVar(dems,xdim); % get the actual data contained
        elseif strcmp(axisAtt,'Y') % if axis attribute is 'Y', the current dimension is the y-dimension
            ydim=j; % define the dimension id
            ydimdata=netcdf.getVar(dems,ydim); % get the actual data contained
        elseif strcmp(axisAtt,'T') % if axis attribute is 'T', the current dimension is the time dimension
            tdim=j; % define the dimension id
            tdimdata=netcdf.getVar(dems,tdim); % get the actual data contained
        end
    end
    % check whether all dimensions exist
    existtest = exist('xdimdata','var');
    existtest = existtest + exist('ydimdata','var');
    existtest = existtest + exist('tdimdata','var');
    % throw errors if not all three dimensions are found
    if existtest==0
        error('All dimesions do not contain axis attribute');
    elseif existtest<3 && existtest>0
        error('Only one or two dimensions contain axis attributes. It should be three');
    end
end

%% Checking how much data is in existing file
if new==false && length(tdimdata)<length(foldersnew) % Data already exists but is not complete
    begin=length(tdimdata)+1;
elseif new==true % No data exists yet
    begin=1;
else % Data is already complete
    netcdf.close(fileid)
    return
end

%% Attributes
% Global attributes
% These are placed every time, so if changes are made, the old attributes
% are overwritten
netcdf.putAtt(fileid,NC_GLOBAL,'title',title) 
netcdf.putAtt(fileid,NC_GLOBAL,'long_title',long_title) 
netcdf.putAtt(fileid,NC_GLOBAL,'comments',comments)
netcdf.putAtt(fileid,NC_GLOBAL,'institution',institution)
netcdf.putAtt(fileid,NC_GLOBAL,'source',source)
netcdf.putAtt(fileid,NC_GLOBAL,'references',references)
netcdf.putAtt(fileid,NC_GLOBAL,'history',history)
netcdf.putAtt(fileid,NC_GLOBAL,'summary',summary)
netcdf.putAtt(fileid,NC_GLOBAL,'keywords',keywords)

netcdf.putAtt(fileid,NC_GLOBAL,'Conventions',conventions)   
netcdf.putAtt(fileid,NC_GLOBAL,'Conventions_help',conventionshelp)

netcdf.putAtt(fileid,NC_GLOBAL,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'))
netcdf.putAtt(fileid,NC_GLOBAL,'CreatedBy',author)
netcdf.putAtt(fileid,NC_GLOBAL,'Correspondence_to',correspondence)

if new==true
    % DEM group attributes
    netcdf.putAtt(dems,NC_GLOBAL,'title',dems_title) 
    netcdf.putAtt(dems,NC_GLOBAL,'long_title',dems_long_title) 
    netcdf.putAtt(dems,NC_GLOBAL,'history',dems_history)
    netcdf.putAtt(dems,NC_GLOBAL,'description',dems_description)

    % Metronome attributes
    % General
    netcdf.putAtt(fileid,NC_GLOBAL,'experiment_nr',expnr)
    netcdf.putAtt(fileid,NC_GLOBAL,'conducted_by',labslaves)
    netcdf.putAtt(fileid,NC_GLOBAL,'metronome_comments',metronome_comments)
    netcdf.putAtt(fileid,NC_GLOBAL,'manmade_structures',manmade_structures)
    netcdf.putAtt(fileid,NC_GLOBAL,'nr_of_cycles',nr_of_cycles)
    netcdf.putAtt(fileid,NC_GLOBAL,'start_of_experiment',start_of_experiment)
    netcdf.putAtt(fileid,NC_GLOBAL,'end_of_experiment',end_of_experiment)
    % Boundary conditions - tilting and weir
    netcdf.putAtt(fileid,NC_GLOBAL,'amplitude_tilt_1',amplitude_tilt_1)
    netcdf.putAtt(fileid,NC_GLOBAL,'period_tilt_1',period_tilt_1)
    netcdf.putAtt(fileid,NC_GLOBAL,'phase_tilt_1',phase_tilt_1)
    netcdf.putAtt(fileid,NC_GLOBAL,'offset_tilt_1',offset_tilt_1)
    netcdf.putAtt(fileid,NC_GLOBAL,'amplitude_tilt_2_overtide',amplitude_tilt_2_overtide)
    netcdf.putAtt(fileid,NC_GLOBAL,'period_tilt_2_overtide',period_tilt_2_overtide)
    netcdf.putAtt(fileid,NC_GLOBAL,'phase_tilt_2_overtide',phase_tilt_2_overtide)
    netcdf.putAtt(fileid,NC_GLOBAL,'amplitude_weir',amplitude_weir)
    netcdf.putAtt(fileid,NC_GLOBAL,'period_weir',period_weir)
    netcdf.putAtt(fileid,NC_GLOBAL,'offset_weir',offset_weir)
    netcdf.putAtt(fileid,NC_GLOBAL,'phase_weir_in_relation_to_tilting',phase_weir_in_relation_to_tilting)
    netcdf.putAtt(fileid,NC_GLOBAL,'amplitude_weir_overtide',amplitude_weir_overtide)
    netcdf.putAtt(fileid,NC_GLOBAL,'period_weir_overtide',period_weir_overtide)
    netcdf.putAtt(fileid,NC_GLOBAL,'phase_weir_overtide',phase_weir_overtide)
    netcdf.putAtt(fileid,NC_GLOBAL,'flood-ebb_dominance',flood_ebb_dominance)
    % Boundary conditions - river
    netcdf.putAtt(fileid,NC_GLOBAL,'river_water_discharge',river_water_discharge)
    netcdf.putAtt(fileid,NC_GLOBAL,'river_water_timing',river_water_timing)
    netcdf.putAtt(fileid,NC_GLOBAL,'river_sand_feed_rate',river_sand_feed_rate)
    netcdf.putAtt(fileid,NC_GLOBAL,'river_sand_timing',river_sand_timing)
    netcdf.putAtt(fileid,NC_GLOBAL,'river_walnut_feed_rate',river_walnut_feed_rate)
    netcdf.putAtt(fileid,NC_GLOBAL,'river_walnut_timing',river_walnut_timing)
    % Boundary conditions - waves
    netcdf.putAtt(fileid,NC_GLOBAL,'waves_on-off',waves_on_off)
    netcdf.putAtt(fileid,NC_GLOBAL,'waves_frequency',waves_frequency)
    netcdf.putAtt(fileid,NC_GLOBAL,'waves_timing',waves_timing)
    % Boundary conditions - sediment at inlet
    netcdf.putAtt(fileid,NC_GLOBAL,'sediment_at_inlet_yes-no',sediment_at_inlet_yes_no)
    netcdf.putAtt(fileid,NC_GLOBAL,'inlet_sand_feed_rate',inlet_sand_feed_rate)
    netcdf.putAtt(fileid,NC_GLOBAL,'inlet_sand_timing',inlet_sand_timing)
    netcdf.putAtt(fileid,NC_GLOBAL,'inlet_walnut_feed_rate',inlet_walnut_feed_rate)
    netcdf.putAtt(fileid,NC_GLOBAL,'inlet_walnut_timing',inlet_walnut_timing)
    % Boundary conditions - sea-level rise
    netcdf.putAtt(fileid,NC_GLOBAL,'sea-level_rise_yes-no',sea_level_rise_yes_no)
    netcdf.putAtt(fileid,NC_GLOBAL,'sea-level_rise_start_cycle',sea_level_rise_start_cycle)
    netcdf.putAtt(fileid,NC_GLOBAL,'sea-level_rise_end_cycle',sea_level_rise_end_cycle)
    netcdf.putAtt(fileid,NC_GLOBAL,'sea-level_rise_cycle_increment',sea_level_rise_cycle_increment)
    netcdf.putAtt(fileid,NC_GLOBAL,'sea-level_rise_weir_change_increment',sea_level_rise_weir_change_increment)
    netcdf.putAtt(fileid,NC_GLOBAL,'sea-level_rise_total_nr_of_events',sea_level_rise_total_nr_of_events)
    % Initial conditions
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_sediment_bed_level',initial_sediment_bed_level)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_coastline_x_position',initial_coastline_x_position)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_channel_depth',initial_channel_depth)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_channel_width',initial_channel_width)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_channel_width_variation',initial_channel_width_variation)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_tidal_basin_depth',initial_tidal_basin_depth)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_tidal_basin_width',initial_tidal_basin_width)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_tidal_basin_length',initial_tidal_basin_length)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_barriers_present',initial_barriers_present)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_barrier_width',initial_barrier_width)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_barrier_height',initial_barrier_height)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_inlet_width',initial_inlet_width)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_barriers_fixated',initial_barriers_fixated)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_landward_barrier_present',initial_landward_barrier_present)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_landward_barrier_x_position',initial_landward_barrier_x_position)
    netcdf.putAtt(fileid,NC_GLOBAL,'initial_landward_barrier_fixated',initial_landward_barrier_fixated)
    % Vegetation
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_yes-no',vegetation_yes_no)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_start_sowing',vegetation_start_sowing)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_interval_sowing',vegetation_interval_sowing)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_nr_of_sowing_events',vegetation_nr_of_sowing_events)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_river_Lotus',vegetation_river_Lotus)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_river_Alfalfa',vegetation_river_Alfalfa)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_river_Veronica',vegetation_river_Veronica)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_inlet_Lotus',vegetation_inlet_Lotus)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_inlet_Alfalfa',vegetation_inlet_Alfalfa)
    netcdf.putAtt(fileid,NC_GLOBAL,'vegetation_inlet_Veronica',vegetation_inlet_Veronica)

else % append history to existing history
    historyfile = netcdf.getAtt(fileid,netcdf.getConstant("NC_GLOBAL"),'history');
    demshistoryfile=netcdf.getAtt(dems,netcdf.getConstant("NC_GLOBAL"),'history');
    historycombined=[historyfile ', ' history];
    demshistorycombined=[demshistoryfile ', ' dems_history];
    netcdf.putAtt(fileid,netcdf.getConstant("NC_GLOBAL"),'history',historycombined)
    netcdf.putAtt(dems,netcdf.getConstant("NC_GLOBAL"),'history',demshistorycombined)
end

%% Processing the laser data
for j=begin:length(foldersnew)
    [DEM(j)] = laserMetronome_ESL([pwd '01metronome_experiments\Exp' exp '\raw_data\laser_scanner\' foldersnew(j).name '\'], expalt, timesteps(j), gridresol, roi, percen);
%     [DEM(j)] = load([pwd '01metronome_experiments\Exp' exp '\processed_data\DEMs\laser_scanner\Exp' expalt '_' foldersnew(j).name '_Grid' num2str(uint8(gridresol*1000)) 'mm']).g;
end

%% Necessary correction outside of laser processing funtion - to be eliminated in the future
corr=-0.0745; % Correction factor valid since 3500 cycles in Exp047
for j=begin:length(foldersnew)
    DEM(j).zlow=DEM(j).zlow+corr;
    DEM(j).z50=DEM(j).z50+corr;
    DEM(j).zhigh=DEM(j).zhigh+corr;
end

%% Optional saving of DEMs as .mat-files
if matsave==true
    for j=begin:length(foldersnew)
        g=DEM(j);
        savename = [pwd '01metronome_experiments\Exp' exp '\processed_data\DEMs\laser_scanner\Exp' expalt '_' foldersnew(j).name '_Grid' num2str(uint8(gridresol*1000)) 'mm'];
        save(savename,'g');
    end
end

%% Adding the data to the NetCDF-file
if new==true % for new files
    % Define the variables and their attributes
    DEMzlow = netcdf.defVar(dems,'zlow','double',[dimidX, dimidY, dimidT]);
    netcdf.putAtt(dems,DEMzlow,'name','zhigh');
    netcdf.putAtt(dems,DEMzlow,'long_name',[num2str(percen(1)) 'th percentile of sediment bed elevation above flume floor']);
    netcdf.putAtt(dems,DEMzlow,'units','m');
    netcdf.putAtt(dems,DEMzlow,'cell_methods',['area: ' num2str(percen(1)) 'th percentile (interval: 1 mm)']);
    netcdf.putAtt(dems,DEMzlow,'valid_range',[0 0.4]);
    netcdf.defVarFill(dems,DEMzlow,false,fillValue);
    
    DEMz50 = netcdf.defVar(dems,'z50','double',[dimidX, dimidY, dimidT]);
    netcdf.putAtt(dems,DEMz50,'name','z50');
    netcdf.putAtt(dems,DEMz50,'long_name','median sediment bed elevation above flume floor');
    netcdf.putAtt(dems,DEMz50,'units','m');
    netcdf.putAtt(dems,DEMz50,'cell_methods','area: median (interval: 1 mm)');
    netcdf.putAtt(dems,DEMz50,'valid_range',[0 0.4]);
    netcdf.defVarFill(dems,DEMz50,false,fillValue);
    
    DEMzhigh = netcdf.defVar(dems,'zhigh','double',[dimidX, dimidY, dimidT]);
    netcdf.putAtt(dems,DEMzhigh,'name','zhigh');
    netcdf.putAtt(dems,DEMzhigh,'long_name',[num2str(percen(2)) 'th percentile of sediment bed elevation above flume floor']);
    netcdf.putAtt(dems,DEMzhigh,'units','m');
    netcdf.putAtt(dems,DEMzhigh,'cell_methods',['area: ' num2str(percen(2)) 'th percentile (interval: 1 mm)']);
    netcdf.putAtt(dems,DEMzhigh,'valid_range',[0 0.4]);
    netcdf.defVarFill(dems,DEMzhigh,false,fillValue);

    DEMn = netcdf.defVar(dems,'n','double',[dimidX, dimidY, dimidT]);
    netcdf.putAtt(dems,DEMn,'name','n');
    netcdf.putAtt(dems,DEMn,'long_name','number of raw data values that was averaged for the grid cell value');
    netcdf.putAtt(dems,DEMn,'units','-');
    netcdf.putAtt(dems,DEMn,'valid_range',[1 1000]);
    netcdf.defVarFill(dems,DEMn,false,fillValue);
else % If file (and variables) already exist
    % Get variable IDs from the file
    DEMzhigh=netcdf.inqVarID(dems,'zhigh');
    DEMzlow=netcdf.inqVarID(dems,'zlow');
    DEMz50=netcdf.inqVarID(dems,'z50');
    DEMn=netcdf.inqVarID(dems,'n');
end

for j=begin:length(foldersnew) % Append data to the variables
    zhighdata=DEM(j).zhigh;
    zhighdata(isnan(zhighdata)) = fillValue;
    zlowdata=DEM(j).zlow;
    zlowdata(isnan(zlowdata)) = fillValue;
    z50data=DEM(j).z50;
    z50data(isnan(z50data)) = fillValue;
    ndata=DEM(j).n;
    ndata(isnan(ndata)) = fillValue;
    netcdf.putVar(dems,DEMzhigh,[0 0 j-1],[nx ny 1],zhighdata'); % Data uses the full x and y dimensions and the length of the time dimension is expanded by one for each DEM
    netcdf.putVar(dems,DEMzlow,[0 0 j-1],[nx ny 1],zlowdata');
    netcdf.putVar(dems,DEMz50,[0 0 j-1],[nx ny 1],z50data');
    netcdf.putVar(dems,DEMn,[0 0 j-1],[nx ny 1],ndata');
    netcdf.putVar(dems,tdim,j-1,1,timesteps(j)); % Adding the cycle number to the time dimension
end
 
%% Closing the file
netcdf.close(fileid)