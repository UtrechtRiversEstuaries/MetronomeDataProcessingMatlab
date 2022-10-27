function [g] = laserMetronome_ESL(rootfolder, Expnr, timestep, gridresol, roi, percen)
% Function to create a Digital Elevation Model on a rectangular grid
% from the line laser scanner data Metronome in Earth Simulation Laboratory
% Maarten G. Kleinhans (December 2022), Utrecht University, Faculty of Geosciences
% 
% produces a structure g (gridded bathymetry data) that contains grids with units in m:
% g.x and g.y: coordinates (m)
% g.z50, g.zlow and g.zhigh: median elevation, low and high percentiles, g.n number of points
%       Note: the Metronome has cartesian coordinates: 
%       x = 0-20 m along the flume to the right
%       y = 0-3 m across the flume to the back wall
%       z = 0-0.4 m elevation above the flume floor
% 
% required input:
% rootfolder: folder where all experiments are stored, e.g. 
%       ['C:\Users\klein101\Documents\Matlabdata\Metronomedata\Exp' Expnr '\']
% Expnr: experiment number, e.g. '046'
% gridresol: required grid resolution (m) (not smaller than scanline spacing 1 mm)
% roi: Region Of Interest as [xmin xmax ymin ymax] (m)
%       if empty then maximised to Metronome dimensions 20 by 3 m
% percen: specify two percentiles e.g. [10 90] for maps, vegetation and outlier removal

%close all
%clear
%clc


%% Calibration settings
%CALIBRATION FOR DECEMBER 2021 Experiment 046 onwards (control sand only for Sarah Hautekiet)
%no calibration in x-direction (controlled by gantry)
scanline            = 0.001;        %m AS SET IN MOTION INTERFACE

%calibration in y-direction for laser-camera system in image pixel coordinates 
opticalcenter       = 1.53;         %m target centered in camera (1.500 is center of the flume)
transversescale     = 0.976*0.001;  %scale number

%calibration in z-direction for laser-camera system in 16 bit image %coordinates within Area of Interest in camera
cameradist          = 5.5;          %strictly: sqrt(4^2+2^2); %4 m height, camera 2 m from laser BUT HERE CALIBRATED
verticalscale       = 1.2987e-4;    %scale number CHECKED ON CALIBRATION TARGET
verticaloffset      = 1498;         %value at flume bottom

%correction of z for curvature (perhaps optical distortion)
parabfactor         =  0.38e-3;     %multiplication
parabhorizoffset    =  1.0;         %m
parabvertioffset    = -0.65e-3;     %m


%% Initialisation
% check input
if nargin == 3
    gridresol = 10/1000; %m
    roi = [0 20 0 3]; %m
    percen = [10 90]; %m
end
if nargin == 4
    roi = [0 20 0 3]; %m
    percen = [10 90]; %m
end
if nargin == 5
    percen = [10 90]; %m
end

%% collect input
files       = dir([rootfolder '*.tif']);
%files       = dir([rootfolder 'Exp' Expnr '_' timestep '*.tif']);
numfiles    = length(files);
testIM      = imread([rootfolder files(1).name]);
IMsize      = size(testIM); %images 4096 wide
clear testIM

%% Coordinates in scan data images
Xi = 0:scanline:20-scanline;
Yi = ( transversescale/2:transversescale:transversescale*IMsize(2) ) ...
    + opticalcenter - transversescale*(IMsize(2))/2;

% Coordinates of grid cell centers in gridded data
Xc = (roi(1):gridresol:roi(2))'; %arrays corner coordinates in x and in y direction
Yc = (roi(3):gridresol:roi(4))';
Xg = (roi(1)+(gridresol/2):gridresol:roi(2)-(gridresol/2))'; %arrays center coordinates in x and in y direction
Yg = (roi(3)+(gridresol/2):gridresol:roi(4)-(gridresol/2))';
[x,y] = meshgrid(Xg,Yg);
g.x = x; g.y = y; %clear x y
g.z50 = NaN(size(g.x)); g.zlow = g.z50; g.zhigh = g.z50; g.n = g.z50; 


%% Load and join laser data from raw tif format
ImDataraw = NaN(IMsize(2),numfiles*IMsize(1));
for i       = 1:numfiles
    %disp(['image ' num2str(i,2) '/' num2str(numfiles,2)]);
    IM      = imread([rootfolder files(i).name]);
    IM      = flipud(rot90(IM,1)); %laser data needs 90 deg rotation and image to real coordinates
    ImDataraw(:,1+((i-1)*IMsize(1)):IMsize(1)+((i-1)*IMsize(1))) = IM;
end
clear IM


%% CALIBRATION
%Calibrate data in vertical
ImDataC = verticalscale.*(ImDataraw-verticaloffset);

%Calibrate coordinates in scan data images
[xi,yiraw] = meshgrid(Xi,Yi);
%Calibrate data in transverse horizontal
yi = (yiraw-opticalcenter) .* (1 - ImDataC./cameradist) +opticalcenter;

%Correct calibrated data for transverse shape
ImData = ImDataC +parabfactor.*(yi -parabhorizoffset).^2 +parabvertioffset;

%Remove unrealistic outliers
ImData( ImData>0.4 | ImData<0 ) = NaN;

%clear ImDataraw ImDataC yiraw


%% Interpolate on grid
waittxt = ['gridding ',num2str(length(Xg)*length(Yg)),' gridcells',' Exp',num2str(Expnr),'-',timestep];
h = waitbar(0,waittxt,'Name','laser Metronome');

for iXg = 1:length(Xg)
    iXi = [ find(Xi>=Xc(iXg),1), find(Xi<=Xc(iXg+1),1,'last') ];
    stripZ = ImData(:,iXi(1):iXi(2)); stripZ = stripZ(:);
    stripY = yi(:,iXi(1):iXi(2));     stripY = stripY(:);
    stripX = xi(:,iXi(1):iXi(2));     stripX = stripX(:);
    for iYg = 1:length(Yg)
        iYi = find( stripY>=Yc(iYg) & stripY<=Yc(iYg+1) );
        tempdat = stripZ(iYi);
        tempdat = tempdat(~isnan(tempdat)); %remove NaN
        ni = length(tempdat);
        g.n(iYg,iXg) = ni;
        if ni>1
            g.z50(iYg,iXg)   = median( tempdat );
            g.zlow(iYg,iXg)  = prctile( tempdat, percen(1));
            g.zhigh(iYg,iXg) = prctile( tempdat, percen(2));
        elseif ni==1
            g.z50(iYg,iXg) = tempdat(1);
        else
            g.z50(iYg,iXg) = NaN;
        end
    end
    waitbar(iXg/length(Xg),h);
end
close(h);


%% Fill the gaps by interpolation
nodatai = find(isnan(g.z50) == 1);
datai   = find(isnan(g.z50) == 0);
g.z50(nodatai)   = griddata(g.x(datai),g.y(datai),g.z50(datai),  g.x(nodatai),g.y(nodatai),'linear');
g.zlow(nodatai)  = griddata(g.x(datai),g.y(datai),g.zlow(datai), g.x(nodatai),g.y(nodatai),'linear');
g.zhigh(nodatai) = griddata(g.x(datai),g.y(datai),g.zhigh(datai),g.x(nodatai),g.y(nodatai),'linear');
