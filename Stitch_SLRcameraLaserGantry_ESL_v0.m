% Code to create overhead images from the SLR camera on the laser gantry
% of the Metronome in Earth Simulation Laboratory
% Maarten G. Kleinhans (APRIL 2022), Utrecht University, Faculty of Geosciences
% 
% produces corrected stitched image and LAB image
% 
%       Note: the Metronome has cartesian coordinates: 
%       x = 0-20 m along the flume to the right
%       y = 0-3 m across the flume to the back wall
%       z = 0-0.4 m elevation above the flume floor
% 
% required input:
% rootfolder: folder where all experiments are stored, e.g. 
%       ['C:\Users\klein101\Documents\Matlabdata\Metronomedata\Exp' Expnr '\']
% Expnr: experiment number, e.g. '046'
% ncam: number of cameras (in case not all were working)
%
% image horizontal is perpendicular to flume, top is upstream
% TO DO:
% geometric correction better with new targets (click their coordinates)
% crop camera images exactly so neither overlap nor gap
% crop entire image to fit Metronome exactly


clear all %variables
close all
clc

%rootfolder  = ['C:\Users\klein101\Documents\Matlabdata\Metronomedata\Exp' num2str(Expnr) '\'];
%rootfolder  = ['C:\Users\Klein101\Documents\Matlabdata\Metronomedata\Exp047overhead\'];
rootfolder  = ['C:\Users\jan-e\OneDrive - Universiteit Utrecht\Uni\Master Thesis\Shared folder salt marsh experiments\SaltMarshMetronome\01raw_metronome_data\Exp049\raw_data\SLR_camera_laser_gantry\10000_dry\'];
Expnr       = 49;   %experiment number
t           = 10000; %tidal cycle
gantrypos = [9.5:0.5:19.5]'; %as set in the motion controller (length same as files)


%% INVENTORY and sort images for automated looping
% list images from rootfolder AND in subdirectories
files = dir([rootfolder '*.jpg']); %/**/ if also underlying directories needed
numfiles = numel(files);
x.cam = zeros(numfiles,1); x.IMnr = x.cam; x.date = x.cam; x.time = x.cam;
% fill list
% for f = 1:numfiles
%     [~, name, ~] = fileparts(files(f).name);
%     filenamestrings = strsplit(name,'_');    
%     x.date(f) = str2double(filenamestrings{end-4});   
%     x.time(f) = str2double(filenamestrings{end-3});    
%     x.datetime(f) = datenum([filenamestrings{end-4} filenamestrings{end-3}],'yymmdd HHMMSS'); %make string that contains date and time
%     %clear filenamestrings
%     %x.filenamec{f} = strcat(name, '.tif');
% end
% clear name;


%% for INTERNAL camera correction: camera and lens parameters
IMsize      = [4000,6000];
% first GUESS of lens calibration, to be determined
params = struct('fx',5000,'fy',5000,'cx',3000,'cy',2000,'k1',-0.11,'k2',0.11,'k3',0,'p1',0,'p2',0,'sk',0);

cameraParams = cameraIntrinsics([params.fx,params.fy],[params.cx,params.cy],IMsize,...
    'RadialDistortion',[params.k1,params.k2,params.k3],...
    'TangentialDistortion',[params.p1,params.p2],'Skew',params.sk);


%% for EXTERNAL camera correction: tform points for projection of camera position
% % positions measured on raw image at 10m at wooden board on the axis of the flume
% xhall =  359; xwall = 5218; 
% yhall = 1781; ywall = 1791; %showing rotation and shift

% positions measured on undistorted image at 10m at 10 m on the tape
xhall2 =  298; xwall2 = 5255;   %showing pixel scale: should be 3000 mm exactly
yhall2 = 2279; ywall2 = 2298;   %showing rotation and shift: should be equal
scalex = (xwall2-xhall2)./3000; %divided by flume width

% positions measured on undistorted image at 10.5m at 10 m on the tape
xhall3 =  295; xwall3 = 5252;   %showing pixel scale: should be 3000 mm exactly
yhall3 = 1445; ywall3 = 1476;   %showing rotation and shift: should be equal
scaley = ( ( ywall2-ywall3 + yhall2-yhall3 )/2 )./ ...
    (1000*(gantrypos(2)-gantrypos(1)) );  %divided by gantry move distance
%scale = (scalex+scaley)/2;
pixshift = round(scaley*(1000*(gantrypos(2)-gantrypos(1))/2 ) ); %gantry movement in pixels

% indices which image parts to keep (REDO AFTER geometric correction)
iwidth  = [290, 5260]; %[295, 5255];
iheight = [IMsize(1)/2-pixshift, IMsize(1)/2+pixshift];
nwidth  = iwidth(2)-iwidth(1)+1;
nheight = iheight(2)-iheight(1)+1;

% DO for one image
% fixedPoints(:,:)  = [xboth(1) yhall; xboth(1) ywall; xboth(2) yhall; xboth(2) ywall];
% movingPoints(:,:) = [819 1815; 828 245; 1947 1836; 1993 267]; %6 16 5 15
% 
% tform = fitgeotrans(squeeze(movingPoints(:,:)),squeeze(fixedPoints(:,:)),'projective');


%% LOOP through corrections
ims(:,:,1) = repmat(uint8(60), nwidth,  numfiles*nheight);
ims(:,:,2) = ims(:,:,1);
ims(:,:,3) = ims(:,:,1);
%figure, imshow(ims)

for p=1:numfiles; %position counter

    % read image
    im = imread([files(p).folder '\' files(p).name]);
    %figure, imshow(im)

    % correct camera distortion (lens etc)
    im2 = undistortImage(im,cameraParams);
    %figure, imshow(im2)

    % re-project
    %im3 = imwarp(im3,tform{c});

    % rotate and select relevant strip from image
    im3{p} = rot90( im2(iheight(1):iheight(2),iwidth(1):iwidth(2),:) ); %also do flip for cartesian coordinates
    %im3 = im2(iheight(1):iheight(2),iwidth(1):iwidth(2),:); %also do flip for cartesian coordinates
    %figure, imshow(im3)

end

for p=1:numfiles; %position counter
    % put into total image
    ims(:, 1+(p-1)*nheight:p*nheight ,:) = im3{p};
end

figure, imshow(ims)
% DONE! This is the image
imwrite(ims, ['C:\Users\jan-e\OneDrive - Universiteit Utrecht\Uni\Master Thesis\Shared folder salt marsh experiments\SaltMarshMetronome\01raw_metronome_data\Exp' num2str(Expnr,'%03i') '\processed_data\orthophotos\SLR_camera_laser_gantry\SLRgantry_Exp' num2str(Expnr,'%03i') '_' num2str(t,'%05i') '_colour.tif']);


%% IMAGE calculations EXAMPLES
%LAB colour space for water depth and vegetation: 
% Luminosity, +Red to -Green, +Yellow to -Blue
% Luminosity is nice grayscale image
% -Green could work to get vegetation, -Blue for water depth
imsLAB = rgb2lab(ims);
figure, hold on, 
subplot(3,1,1), imshow(imsLAB(:,:,1),[0 70]), colorbar, title('luminosity')
subplot(3,1,2), imshow(imsLAB(:,:,2),[-30 -10]); colorbar, title('green-ishness')
subplot(3,1,3), imshow(imsLAB(:,:,3),[-10 20]); colorbar, title('blue-ishness')
% print('LABchannels','-dpng');

% convert back RGB with enhanced colours
imsLABc = imsLAB;
colourenhance = 1.5; %whatever you like but don't exceed max range
imsLABc(:,:,2) = colourenhance.*imsLAB(:,:,2); %multiply the colour depth
imsLABc(:,:,3) = colourenhance.*imsLAB(:,:,3);
imsRGBc = lab2rgb(imsLABc); 
figure, imshow(imsRGBc); %coloured like a tropical system :-)

% DONE! This is the image; replace colour by blueness etc
imwrite(imsRGBc, ['C:\Users\jan-e\OneDrive - Universiteit Utrecht\Uni\Master Thesis\Shared folder salt marsh experiments\SaltMarshMetronome\01raw_metronome_data\Exp' num2str(Expnr,'%03i') '\processed_data\orthophotos\SLR_camera_laser_gantry\SLRgantry_Exp' num2str(Expnr,'%03i') '_' num2str(t,'%05i') '_colourenhanced.tif']);
% if GeoTIFF needed, get coordinates right after cropping:
% ims = flip(rot90(ims));
% R = georefcells([0 20],[0 3],size(ims));
% geotiffwrite(['STITCH\Exp' num2str(Expnr,'%03i') '_RGB_' num2str(t,'%04i') '.tif'],foto1,R);

