% Code to create overhead images and movies
% from the overhead cameras of the Metronome in Earth Simulation Laboratory
% Maarten G. Kleinhans (APRIL 2022), Utrecht University, Faculty of Geosciences
% 
% produces corrected stitched image and blueness image, can produce movie
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
% TO DO:
% geometric correction better with new targets (click their coordinates)
% crop camera images exactly so neither overlap nor gap
% colour corrections
% crop entire image to fit Metronome exactly


clear all %variables
close all
clc

Expnr       = 49;
%rootfolder  = ['C:\Users\klein101\Documents\Matlabdata\Metronomedata\Exp' num2str(Expnr) '\'];
%rootfolder  = ['C:\Users\Klein101\Documents\Matlabdata\Metronomedata\Exp047overhead\'];
rootfolder  = ['Z:\Metronome\experiments\Exp049\raw_data\overhead_cameras'];
ncam = 4; %number of cameras


%% INVENTORY and sort images for automated looping
% list images from rootfolder AND in subdirectories
files = dir([rootfolder '/**/*.bmp']);
numfiles = numel(files);
% if numfiles==0
%     error('no images, or not located in timelapse folder');
% end
% prepare list
x.cam = zeros(numfiles,1); x.IMnr = x.cam; x.date = x.cam; x.time = x.cam;
% fill list
for f = 1:numfiles
    [~, name, ~] = fileparts(files(f).name);
    filenamestrings = strsplit(name,'_');    
    x.cam(f)  = str2double(filenamestrings{1});    
    x.IMnr(f) = str2double(filenamestrings{2});    
    x.date(f) = str2double(filenamestrings{3});   
    x.time(f) = str2double(filenamestrings{4});    
    x.datetime(f) = datenum([filenamestrings{3} filenamestrings{4}],'yyyymmdd HHMMSSFFF'); %make string that contains date and time
    clear filenamestrings
    x.filenamec{f} = strcat(name, '.tif');     
end
clear name;

% sort by camera and date
nr      = zeros(ncam,1);   %allocate
anyNaN  = cell(ncam,1);
% icam    = zeros(round(1.1*numfiles/ncam),ncam); %allocate for estimated timesteps
icam    = zeros(round(numfiles/ncam),ncam); %allocate for estimated timesteps
time    = icam;
for c   = 1:ncam
    j       = find(x.cam==c);
    nr(c)   = numel(j);  
    NR      = sum(nr);
    temp(:,1) = j(:);
    temp(:,2) = x.datetime(NR-nr(c)+1:NR);
    out     = sortrows(temp,2);
    %j = find(out(:,2)==0);
    %out(j,:) = [];
    icam(1:length(out),c) = out(:,1); %needs at least two timesteps
    time(1:length(out),c) = out(:,2);
    clear temp out j
end

% check visually for missing images
% plot(time,'+')

%% for INTERNAL camera correction: camera and lens parameters
IMsize      = [2048,2048]; %images 2048 by 2048
% params(1) = struct('fx',2305.83068350631,'fy',2305.83068350631,'cx',1040.49151827290,'cy',1025.97769732717,'k1',-0.200057651376996,'k2',0.180513910670662,'k3',-0.0587050208963713,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.534174540670714);
% params(2) = struct('fx',2311.49278149317,'fy',2311.49278149317,'cx',1004.65655498067,'cy',1002.72726837942,'k1',-0.202745988954746,'k2',0.170838989656136,'k3', 0.0331307105717511,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.470337051946155);
% params(3) = struct('fx',2317.46492335671,'fy',2317.46492335671,'cx',1029.18847346011,'cy',1004.86012699160,'k1',-0.203881376570696,'k2',0.208025873140558,'k3',-0.0912079096853060,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.309785152058378);
% params(4) = struct('fx',2305.46599866940,'fy',2305.46599866940,'cx',1013.11645024591,'cy',1028.73327343962,'k1',-0.194492771474653,'k2',0.130376993700372,'k3', 0.0968288910295995,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.486141852568340);
% params(5) = struct('fx',2309.38451783626,'fy',2309.38451783626,'cx',1019.59383023124,'cy',1027.36476114956,'k1',-0.199737456248573,'k2',0.161255010924189,'k3', 0.0270063891072896,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.280646555345602);
% params(6) = struct('fx',2315.41010096867,'fy',2315.41010096867,'cx',1020.29127418637,'cy',1021.36319641566,'k1',-0.198979859453963,'k2',0.188914047408094,'k3',-0.1446060809629000,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.459454526216367);
% params(7) = struct('fx',2312.98557382207,'fy',2312.98557382207,'cx',1006.73388170178,'cy',1033.58969878715,'k1',-0.195789963523382,'k2',0.124103122428344,'k3', 0.1522651103649300,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.554492875721716);
% params(8) = struct('fx',2311,'fy',2311,'cx',1024,'cy',1024,'k1',-0.2,'k2',0.17,'k3',0.1,'p1',0,'p2',0,'sk',0);

% FOUR cameras 4 to 7 for Exp47 after grabber card failure
params(1) = struct('fx',2305.46599866940,'fy',2305.46599866940,'cx',1013.11645024591,'cy',1028.73327343962,'k1',-0.194492771474653,'k2',0.130376993700372,'k3', 0.0968288910295995,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.486141852568340);
params(2) = struct('fx',2309.38451783626,'fy',2309.38451783626,'cx',1019.59383023124,'cy',1027.36476114956,'k1',-0.199737456248573,'k2',0.161255010924189,'k3', 0.0270063891072896,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.280646555345602);
params(3) = struct('fx',2315.41010096867,'fy',2315.41010096867,'cx',1020.29127418637,'cy',1021.36319641566,'k1',-0.198979859453963,'k2',0.188914047408094,'k3',-0.1446060809629000,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.459454526216367);
params(4) = struct('fx',2312.98557382207,'fy',2312.98557382207,'cx',1006.73388170178,'cy',1033.58969878715,'k1',-0.195789963523382,'k2',0.124103122428344,'k3', 0.1522651103649300,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.554492875721716);

for c=1:ncam
    cameraParams(c) = cameraIntrinsics([params(c).fx,params(c).fy],[params(c).cx,params(c).cy],IMsize,...
        'RadialDistortion',[params(c).k1,params(c).k2,params(c).k3],...
        'TangentialDistortion',[params(c).p1,params(c).p2],'Skew',params(c).sk);
end

% illumination tests date 20211014
illumfac = 1+[-0.049353916	0.041177635	-0.010513378	0.039183768	-0.006236532	-0.010214298	-0.00404327];
% colour balance tests date 20211008 Red/Green and Blue/Green
colourfac = 1-[-0.018803419	-0.01804758;
    -0.006837607	-0.006562756;
    -0.006837607	-0.006562756;
    0.011111111	0.022149303;
    0.005128205	-0.012305168;
    0.011111111	0.027891715;
    0.005128205	-0.006562756 ];


%% for EXTERNAL camera correction: tform points for projection of camera position
% control points wallside: 6 5 4 3 7 2 1, hallside: 16 15 14 13 17 12 11, use millimeters
yhall = 160; ywall = 2850;
xboth = [9.80 11.80 13.50 15.80 16.80 17.50 19.80].*1000;

fixedPoints(1,:,:)  = [xboth(1) yhall; xboth(1) ywall; xboth(2) yhall; xboth(2) ywall];
fixedPoints(2,:,:)  = [xboth(2) yhall; xboth(2) ywall; xboth(3) yhall; xboth(3) ywall];
fixedPoints(3,:,:)  = [xboth(4) yhall; xboth(4) ywall; xboth(5) yhall; xboth(5) ywall];
fixedPoints(4,:,:)  = [xboth(6) yhall; xboth(6) ywall; xboth(7) yhall; xboth(7) ywall];
movingPoints(1,:,:) = [819 1815; 828 245; 1947 1836; 1993 267]; %6 16 5 15
movingPoints(2,:,:) = [335 1824; 326 261; 1295 1826; 1309 268]; %5 15 4 14
movingPoints(3,:,:) = [979 1828; 990 262; 1554 1820; 1547 267]; %3 13 7 17
movingPoints(4,:,:) = [476 1841; 470 238; 1609 1790; 1646 286]; %2 12 1 11

tform = cell(1,ncam); tform{1,ncam} = []; %allocate
for c = 1:ncam
    tform{c} = fitgeotrans(squeeze(movingPoints(c,:,:)),squeeze(fixedPoints(c,:,:)),'projective');
end


%% LOOP through corrections
skip = 250; %number of timesteps to skip

for t=1:skip:nr(1) %time loop

    for c=1:ncam
        % de-Bayer to obtain colour image
        im = demosaic(imread([files(icam(t,c)).folder '\' files(icam(t,c)).name]),'gbrg');

        % correct camera distortion (lens etc)
        im3 = undistortImage(im,cameraParams(c)); %im2

        % re-project
        im4 = imwarp(im3,tform{c});

        % put images together
        im5{c} = im3; %im4 to include warp correction

        % correct illumination from lab tests, FOUR cameras so correct the counting
        im5{c} = illumfac(c+3).*im5{c};

        % correct colour with white balance from lab tests, relative to Green, FOUR cameras so correct the counting
        im5{c}(:,:,1) = colourfac(c+3,1).*im5{c}(:,:,1); %Red band
        im5{c}(:,:,3) = colourfac(c+3,2).*im5{c}(:,:,3); %Blue band

        % % Further manual corrections to tweak to perfection as follows:
        % im5{camera}(:,:,band) = imadjust(im5{camera}(:,:,band),[0 0.9],[0 1]);

    end %of camera loop

    % STITCH cameras into one large image
    % manual positioning, TO IMPROVE X-range after better geometry, crop to fit exactly to flume
    ims = [squeeze(im5{1,1}(96:1995,601:end-200,:)) ...
           squeeze(im5{1,2}(96:1995,201:end-200,:)) ...
           squeeze(im5{1,3}(96:1995,201:end-200,:)) ...
           squeeze(im5{1,4}(96:1995,201:end-200,:)) ...
           ]; %without warp correction
    
    % DONE! This is the image; replace colour by blueness etc
    imwrite(ims, [rootfolder '\overhead_Exp' num2str(Expnr,'%03i') '_' num2str(t,'%05i') '_colour.tif']);
    
    % if GeoTIFF needed, get coordinates right after cropping:
%     ims = flip(rot90(ims));
%     R = georefcells([0 20],[0 3],size(ims));
%     geotiffwrite(['STITCH\Exp' num2str(Expnr,'%03i') '_RGB_' num2str(t,'%04i') '.tif'],foto1,R);


end %of time loop


%% plot examples
%figure, imshow(im3)
%figure(1), imshow(im5{1}), figure(2), imshow(im5{2}), figure(3), imshow(im5{3}), figure(4), imshow(im5{4})
%figure, imshow(im2), figure, imshow(imC)
figure, imshow(ims); 


%% IMAGE calculations EXAMPLES
% turn into doubles and plot
% %imdat = double(ims); 
% %figure, surf(imdat(:,:,1)), shading interp, set(gca,'dataaspectratio',[1 1 1],'view',[0,90]), colorbar

%LAB colour space for water depth and vegetation: 
% Luminosity, +Red to -Green, +Yellow to -Blue
% Luminosity is nice grayscale image
% -Green could work to get vegetation, -Blue for water depth
imsLAB = rgb2lab(ims);

% %get the range of Digital Numbers in the LAB bands
% temp = imsLAB(:,:,1); rangeL = [min(temp(:)) max(temp(:))];
% temp = imsLAB(:,:,2); rangeA = [min(temp(:)) max(temp(:))];
% temp = imsLAB(:,:,3); rangeB = [min(temp(:)) max(temp(:))];

figure, hold on, 
subplot(3,1,1), imshow(imsLAB(:,:,1),[0 70]), colorbar, title('luminosity')
subplot(3,1,2), imshow(imsLAB(:,:,2),[-30 -10]); colorbar, title('green-ishness')
subplot(3,1,3), imshow(imsLAB(:,:,3),[-10 20]); colorbar, title('blue-ishness')

% % remove reflection of darker central laser opening in the illumination
% % TO IMPROVE this on the entire initial, featureless bed (here only clean sand bed chosen)
% lumos = mean(imsLAB(:,1030:1330,1)')'; plot(lumos) 
% %about 873 to 916 is a parabola shaped shadow 10% darker at center 896
% parabx = (-21:21)'; parabfac = 1+0.05.*(1-(parabx./21).^2); %plot(parabx,parabfac)
% imsLABs = imsLAB; 
% imsLABs(896-21:896+21,:,1) = parabfac.*imsLAB(896-21:896+21,:,1);
% figure, imshow(imsLABs(:,:,1),rangeL); colorbar %luminosity only

% convert back RGB with enhanced colours
imsLABc = imsLAB;
colourenhance = 1.5; %whatever you like but don't exceed max range
imsLABc(:,:,2) = colourenhance.*imsLAB(:,:,2); %multiply the colour depth
imsLABc(:,:,3) = colourenhance.*imsLAB(:,:,3);
imsRGBc = lab2rgb(imsLABc); 
figure, imshow(imsRGBc); %coloured like a tropical system :-)


%% MAARTEN's CODE UNTIL HERE











%% OLD STUFF Lisanne Steven, use to plunder code
%% Switches for choices
tic

movieon_rgb = 0;        % RGB movie
movieon_lab = 0;        % Blue movie
rgb_img     = 0;        % RGB picture tif
rgb_fig     = 1;        % RGB png (smaller)
lab_img     = 0;        % blue tif in grayscale
lab_fig     = 0;        % blue png in grayscale

Expnr       = 20211105;
waln        = 0;        % 1 = yes nutshell, 0 = no
int_mov     = 50;
int_img     = 50;       % choice, or same as int_mov
int_fig     = 1;

bestcam     = 6;        % camera with the least time gaps

addrim   = 500;   % offset??


%% Make directories in experiment folder
cd(['C:\Users\klein101\Documents\Matlabdata\Metronomedata\Exp' num2str(Expnr,'%03i')]);
addpath 'R:\Metronome\experiments\matlab_lisanne'


if ~exist(['stitch' num2str(Expnr) '_2.mat'],'file')
    mkdir STITCH
    mkdir MOVIE
    mkdir FIGS
    
    %% load files other scripts
    %load('Z:\Metronome\matlab\overheadcams_matlab_lisanne\Colormaps\mymap.mat','mymap')
    load('Z:\Metronome\matlab\overheadcams_matlab_lisanne\Vignetting_correction\scaleMatVignetting.mat','scaleMatVignetting_gray');       
    load('Z:\Metronome\matlab\overheadcams_matlab_lisanne\Color_calibration\color_corr.mat','meanGray','meanR','meanG','meanB');
    
    %% camera parameters
    params(8) = struct('fx',2311,'fy',2311,'cx',1024,'cy',1024,'k1',-0.2,'k2',0.17,'k3',0.1,'p1',0,'p2',0,'sk',0);
    params(1) = struct('fx',2305.83068350631,'fy',2305.83068350631,'cx',1040.4915182729,'cy',1025.97769732717,'k1',-0.200057651376996,'k2',0.180513910670662,'k3',-0.0587050208963713,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.534174540670714);
    params(2) = struct('fx',2311.49278149317,'fy',2311.49278149317,'cx',1004.65655498067,'cy',1002.72726837942,'k1',-0.202745988954746,'k2',0.170838989656136,'k3',0.0331307105717511,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.470337051946155);
    params(3) = struct('fx',2317.46492335671,'fy',2317.46492335671,'cx',1029.18847346011,'cy',1004.8601269916,'k1',-0.203881376570696,'k2',0.208025873140558,'k3',-0.091207909685306,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.309785152058378);
    params(4) = struct('fx',2305.4659986694,'fy',2305.4659986694,'cx',1013.11645024591,'cy',1028.73327343962,'k1',-0.194492771474653,'k2',0.130376993700372,'k3',0.0968288910295995,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.48614185256834);
    params(5) = struct('fx',2309.38451783626,'fy',2309.38451783626,'cx',1019.59383023124,'cy',1027.36476114956,'k1',-0.199737456248573,'k2',0.161255010924189,'k3',0.0270063891072896,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.280646555345602);
    params(6) = struct('fx',2315.41010096867,'fy',2315.41010096867,'cx',1020.29127418637,'cy',1021.36319641566,'k1',-0.198979859453963,'k2',0.188914047408094,'k3',-0.1446060809629,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.459454526216367);
    params(7) = struct('fx',2312.98557382207,'fy',2312.98557382207,'cx',1006.73388170178,'cy',1033.58969878715,'k1',-0.195789963523382,'k2',0.124103122428344,'k3',0.152265110364930,'p1',0.0000860004594832993,'p2',0.0000995314003703556,'sk',-0.554492875721716);
    for c = 1:7
        params(c).cx = params(c).cx+addrim;
        params(c).cy = params(c).cy+addrim;
    end

    %% geometric transformation - aanpassen als cameras worden verhangen
fixedPoints(1,:,:)  = [1000 0; 1000 3000; 3000 0; 3000 3000];
movingPoints(1,:,:) = [473 1921; 2513 1917; 480 551; 2516 564];
fixedPoints(2,:,:)  = [3000 0; 3000 3000; 5000 0; 5000 3000];
movingPoints(2,:,:) = [479 2321; 2508 2354; 492 958; 2540 989];
fixedPoints(3,:,:)  = [6000 0; 6000 3000; 8000 0; 8000 3000];
if Expnr>10
    movingPoints(3,:,:) = [469 2266; 2456 2221; 440 915; 2458 899];
else
    movingPoints(3,:,:) = [499 2261; 2490 2233; 475 902; 2497 895];
end
fixedPoints(4,:,:)  = [9000 0; 9000 3000; 11000 0; 11000 3000];
movingPoints(4,:,:) = [522 2269; 2559 2212; 497 859; 2562 873];
fixedPoints(5,:,:)  = [12000 0; 12000 3000; 14000 0; 14000 3000];
movingPoints(5,:,:) = [471 2156; 2573 2202; 483 803; 2563 741];
fixedPoints(6,:,:)  = [15000 0; 14000 3000; 16000 0; 16000 3000];
movingPoints(6,:,:) = [555 1976; 2535 2573; 527 1308; 2540 1268];
fixedPoints(7,:,:)  = [18000 0; 18000 3000; 20000 0; 20000 3000];
movingPoints(7,:,:) = [491 1861; 2539 1860; 497 484; 2546 508];

    tform = cell(1,7); %allocate
    tform{1,7} = []; %allocate
    for c = 1:7
        tform{c} = fitgeotrans(squeeze(movingPoints(c,:,:)),squeeze(fixedPoints(c,:,:)),'projective');
    end
 
    %% filter names ed.
    files = dir('timelapse_steven\*.bmp');
    numfiles = numel(files);
    if numfiles==0
        error('no images, or not located in timelapse folder');
    end
    x.cam = zeros(numfiles,1); %allocate
    x.IMnr = x.cam; x.date = x.cam; x.time = x.cam; %allocate
    for f = 1:numfiles
        [~, name, ~] = fileparts(files(f).name);
        filenamestrings = strsplit(name,'_');    
%         if filenamestrings{1} == '1' && str2double(filenamestrings{2}) > 11937
%             %skip incorrectly named files
%         elseif filenamestrings{1} == '7' && str2double(filenamestrings{2}) > 8151
%             %skip incorrectly named files
%         else
        x.cam(f) = str2double(filenamestrings{1});    
        x.IMnr(f) = str2double(filenamestrings{2});    
        x.date(f) = str2double(filenamestrings{3});   
%         x.date(f) = datenum([filenamestrings{3} filenamestrings{4}]);
        x.datetime(f) = datenum([filenamestrings{3} filenamestrings{4}],'yyyymmdd HHMMSSFFF'); %make string that contains date and time
        x.time(f) = str2double(filenamestrings{4});    
        clear filenamestrings
        x.filenamec{f} = strcat(name, '.tif');     
%         end
    end
    clear name;

    %%
    nr      = zeros(7,1);   %allocate
    anyNaN  = cell(7,1);
    icam    = zeros(2e4,7); %allocate
    time    = zeros(2e4,7);
    for c   = 1:7
        j       = find(x.cam==c);
        nr(c)   = numel(j);  
        NR      = sum(nr);
        
        % sort photos by date
        temp(:,1) = j(:);
        temp(:,2) = x.datetime(NR-nr(c)+1:NR);
        out     = sortrows(temp,2);
        j = find(out(:,2)==0);
        out(j,:) = [];
        icam(1:length(out),c) = out(:,1);
        time(1:length(out),c) = out(:,2);
        clear temp out j
    end
    
    % Correct for gaps in recorded photos -- assume bestcam has the most data
    % --> bestcam has to go first in the for loop
    order   = 1:7;
    order(order==bestcam) = [];
    order   = [bestcam order];
    for c   = order
        if c    == bestcam
            j = find(icam(:,bestcam)==0); %allocated matrix is larger than data, so get rid of excess rows of zeros
            icam(j,:) = [];
            time(j,:) = [];
            
            %manual time corrections - found in a different script
%             time_errorz = [101 107 126 196 1734 1749 3500-24 5000-34 9000-44];
%             time_gapsize = [5 5 2 6 6 2 10 10 10];
            
%             clear icam time
%             icam(1:9960,1) = 1;
%             time = icam;
% exp039
%             time_errorz = [1 5001-10 7001-20 9001-30];
%             time_gapsize = [10 10 10 10];
% exp023            
            time_errorz = 3001;
            time_gapsize = 100;
            
            for i = 1:length(time_errorz)
                icam(time_errorz(end-(i-1))+time_gapsize(end-(i-1)):end+time_gapsize(end-(i-1)),bestcam) = icam(time_errorz(end-(i-1)):end,bestcam);
                icam(time_errorz(end-(i-1)):time_errorz(end-(i-1))+time_gapsize(end-(i-1))-1,bestcam) = NaN;
                time(time_errorz(end-(i-1))+time_gapsize(end-(i-1)):end+time_gapsize(end-(i-1)),bestcam) = time(time_errorz(end-(i-1)):end,bestcam);
                time(time_errorz(end-(i-1)):time_errorz(end-(i-1))+time_gapsize(end-(i-1))-1,bestcam) = NaN;
            end
            %Now, ref camera is fixed. Now, match the dates of photos of other cameras to those of ref camera.
        else
        
            tel = nr(c);
            for i = 1:length(time)
                if isnan(time(i,bestcam)) %if gap in reference data, also gap in to-be-corrected data
                    tel = tel+1;
                    time(i+1:tel,c) = time(i:tel-1,c);
                    time(i,c) = NaN;
                    icam(i+1:tel,c) = icam(i:tel-1,c);
                    icam(i,c) = NaN;
                elseif ~any(time(i,c))
                    time(i,c) = NaN;
                    icam(i,c) = NaN;
                else
                    if abs(time(i,bestcam)-time(i,c))> 2/(24*3600) %if time-difference with NoGap array is larger than 2 seconds, then photo was not taken.
                    tel = tel+1;
                    time(i+1:tel,c) = time(i:tel-1,c);        
                    time(i,c) = NaN;
                    icam(i+1:tel,c) = icam(i:tel-1,c);        
                    icam(i,c) = NaN;
                    end
                end
            end
        end
        anyNaN{c} = find(isnan(icam(:,c)));
    end
    %%
    stop = min(nr);

    save(['stitch' num2str(Expnr) '_2'],'files','icam','mymap','nr','params','stop','t_concord','x','scaleMatVignetting_gray','meanGray','meanR','meanG','meanB','time','anyNaN');
else
    load(['stitch' num2str(Expnr) '_2']);
end
 
%% Movie
%int = 1;
if movieon_rgb == 1
    if waln == 1
        moviename = ['R:\Metronome\experiments\Exp' num2str(Expnr,'%03i') '_waln\MOVIE\Exp' num2str(Expnr,'%03i') '_RGB_int' num2str(int_mov) '.mp4'];
    else
        moviename = ['R:\Metronome\experiments\Exp' num2str(Expnr,'%03i') '\MOVIE\Exp' num2str(Expnr,'%03i') '_RGB_int' num2str(int_mov) '.mp4'];
    end
  
    vid = VideoWriter(moviename, 'Motion JPEG AVI');  %duurt lang
    vidObj.Quality = 100; % in %
    vid.FrameRate = 6; % in /s
    open(vid);
end

if movieon_lab == 1
    if waln == 1
        moviename2 = ['R:\Metronome\experiments\Exp' num2str(Expnr,'%03i') '_waln\MOVIE\Exp' num2str(Expnr,'%03i') '_LAB_int' num2str(int_mov) '.mp4'];
    else
        moviename2 = ['R:\Metronome\experiments\Exp' num2str(Expnr,'%03i') '\MOVIE\Exp' num2str(Expnr,'%03i') '_LAB_int' num2str(int_mov) '.mp4'];
    end
    vid2 = VideoWriter(moviename2, 'Motion JPEG AVI');  %duurt lang
    vidObj.Quality = 100; % in %
    vid2.FrameRate = 6; % in /s
    open(vid2);
end


%%
countfig_rgb = 0;
countmov_rgb = 0;
countfig_lab = 0;
countmov_lab = 0;
timelapsename = 'timelapse_steven\';
im2 = zeros(2048,2048,3); %allocate
for t = 5:int_img:length(icam)
    im7 = cell(1,7); %allocate
    for c = 1:7
        im6 = zeros([2048+2*addrim 2048+2*addrim 3]); %allocate
        
        %if all photos are present (so, no gaps), do ordinary combination
        if ~any(find(anyNaN{c}(:)==t:-1:t-4)) 
%         im = demosaic(imread(['timelapse_BI_03500\' files(icam{c}(t)).name]),'gbrg');
                im = imlincomb(...
                1/5, demosaic(imread([timelapsename files(icam(t,c)).name]),'gbrg'),...
                1/5, demosaic(imread([timelapsename files(icam(t-1,c)).name]),'gbrg'),...
                1/5, demosaic(imread([timelapsename files(icam(t-2,c)).name]),'gbrg'),...
                1/5, demosaic(imread([timelapsename files(icam(t-3,c)).name]),'gbrg'),...
                1/5, demosaic(imread([timelapsename files(icam(t-4,c)).name]),'gbrg') );
%                 figure;imshow(im);
                
        %if at least 2 photos are present, do different combination
        else
            [~,temp] = find(anyNaN{c}(:)==t:-1:t-4);
            temp2 = 1:5;
            temp2(temp) = [];
            if (numel(temp2)<5 && numel(temp2)>1)
                p = numel(temp2);
                
                commandStr = 'im = imlincomb(';
                for i = 1:length(temp2)
                    if i ~= p
                        commandStr  = strcat(commandStr,['1/' num2str(p) ', demosaic(imread([timelapsename files(icam(t-' num2str(temp2(i)-1) ',c)).name]),''gbrg''),']); 
                    else
                        commandStr  = strcat(commandStr,['1/' num2str(p) ', demosaic(imread([timelapsename files(icam(t-' num2str(temp2(i)-1) ',c)).name]),''gbrg'')']); 
                    end
                end
                commandStr = strcat(commandStr,');');
                eval(commandStr);
%                 figure;imshow(im);
            
            else %1 or all images are missing
                ref     = t-2;
                up      = ref;
                dwn     = ref;
                
                while isnan(icam(up,c))
                    up  = up-1;
                    if up == 0
                        break
                    end
                end
                while isnan(icam(dwn,c))
                    dwn = dwn+1;
                end
                
                if up > 1
                    GapSize = dwn-(up+4); %+4, because smaller gaps are handled by first part of this if-else statement

                    im = imlincomb(...
                    (2-2*(t-up-4)/GapSize)/4, demosaic(imread([timelapsename files(icam(up,c)).name]),'gbrg'),...
                    (2-2*(t-up-4)/GapSize)/4, demosaic(imread([timelapsename files(icam(up-1,c)).name]),'gbrg'),...
                    (0+2*(t-up-4)/GapSize)/4, demosaic(imread([timelapsename files(icam(dwn,c)).name]),'gbrg'),...
                    (0+2*(t-up-4)/GapSize)/4, demosaic(imread([timelapsename files(icam(dwn+1,c)).name]),'gbrg') );
    %                 clf; imshow(im); drawnow; pause(0.1); clc;
    %                 pause(0.5)
    %                 end
                elseif up == 0 %in case the first image(s) are missing, take first nonNaN images
                    while isnan(icam(dwn,c))
                        dwn = dwn+1;
                    end
                    im = imlincomb(...
                    1/2, demosaic(imread([timelapsename files(icam(dwn,c)).name]),'gbrg'),...
                    1/2, demosaic(imread([timelapsename files(icam(dwn+1,c)).name]),'gbrg') );
%                 elseif dwn == 0 %in case the last image(s) are missing, take last nonNaN images
%                     while isnan(icam(dwn,c))
%                         dwn = dwn+1;
%                     end
%                     im = imlincomb(...
%                     1/2, demosaic(imread(['timelapse\' files(icam(dwn,c)).name]),'gbrg'),...
%                     1/2, demosaic(imread(['timelapse\' files(icam(dwn+1,c)).name]),'gbrg') );
                end
            end
        end
        %% remove striping
        % I do not observe lines in cam1,3,5,6,7, a few lines in cam2, and
        % a lot of line in cam4.
        win = 15;
        lineI = zeros(2048,3); bgI = lineI;
        for tel=1:3
            if c==1
                lineI(:,tel) = mean(im(1:1200,:,tel)); %lineavg2 = mean(imavg');
            else
                lineI(:,tel) = mean(im(:,:,tel)); %lineavg2 = mean(imavg');
            end
            bgI(:,tel) = smooth(lineI(:,tel),win,'rlowess')'; %background rlowess or rloess
        end
        resI = lineI - bgI;
        multI = 0.*im;
        for tel=1:3
            multI(:,:,tel) = repmat(resI(:,tel)',2048,1);
        end
        im1 = im - multI;
        % figure;imshow(im1);
        
        %% vignetting
        for k = 1:3
            im2(:,:,k) = double(im1(:,:,k)) .* squeeze(scaleMatVignetting_gray(c,:,:));
        end
        % figure;imshow(uint8(im2));
        
        %% contrast & color
%         im11(:,:,1) = adapthisteq(im1(:,:,1));
%         im11(:,:,2) = adapthisteq(im1(:,:,2));
%         im11(:,:,3) = adapthisteq(im1(:,:,3));
        
%         im3(:,:,1) = imadjust(im2(:,:,1),[0.07 0.35]);
%         im3(:,:,2) = imadjust(im2(:,:,2),[0.07 0.35]);
%         im3(:,:,3) = imadjust(im2(:,:,3),[0.07 0.35]);
        % figure;imshow(im3);
        
        % im3 = imadjust(im2,[0 adj(c,3)]); % this is not working...
        
        % Make all channels have the same mean
        if Expnr<24
            redChannel = uint8(im2(:,:,1) * double(meanGray(c)) / double(meanR(c)));
            greenChannel = uint8(im2(:,:,2) * double(meanGray(c)) / double(meanG(c)));
            blueChannel = uint8(im2(:,:,3) * double(meanGray(c)) / double(meanB(c)));
            % Recombine separate color channels into a single, true color RGB image.
            im3 = cat(3, redChannel, greenChannel, blueChannel);
        else
            im3=uint8(im2);
        end
        % imshowpair(uint8(im2),im3,'montage')
       
        %% rectify
        im4 = im2double(im3);
        im5 = zeros(size(im4,1)+2*addrim,size(im4,2)+2*addrim,size(im4,3));
        im5(addrim+1:addrim+size(im4,1),addrim+1:addrim+size(im4,2),:) = im4;
        % figure;imshow(im5);
        
        for band = 1:3
            im6(:,:,band) = undistort(squeeze(im5(:,:,band)), params(c));
        end
        % figure;imshow(im6);
        
        %% warp
        %cpselect(squeeze(J(4,:,:,2)),squeeze(J(4,:,:,2)))
        im7{c} = imwarp(im6,tform{c});
        % figure;imshow(im6);
        %im5{c} = imtransform(im4,t_concord,'XYScale',[1 1]);
        %figure;subplot(1,2,1);imshow(im4);subplot(1,2,2);imshow(im5{c});
        
    end
    
    im8 = im7;    
    if Expnr<24
        im8{1,1}(:,:,1) = imadjust(double(im8{1,1}(:,:,1)),[0 0.77],[0 1]);
        im8{1,2}(:,:,1) = imadjust(double(im8{1,2}(:,:,1)),[0 0.9],[0 1]);
        im8{1,3}(:,:,1) = imadjust(double(im8{1,3}(:,:,1)),[0 0.875],[0 1]);
        im8{1,4}(:,:,1) = imadjust(double(im8{1,4}(:,:,1)),[0 0.845],[0 1]);
        im8{1,5}(:,:,1) = imadjust(double(im8{1,5}(:,:,1)),[0 0.82],[0 1]);
        im8{1,6}(:,:,1) = imadjust(double(im8{1,6}(:,:,1)),[0 0.75],[0 1]);
        im8{1,7}(:,:,1) = imadjust(double(im8{1,7}(:,:,1)),[0 0.76],[0 1]);   

        im8{1,1}(:,:,2) = imadjust(double(im8{1,1}(:,:,2)),[0 0.78],[0 1]);
        im8{1,2}(:,:,2) = imadjust(double(im8{1,2}(:,:,2)),[0 0.88],[0 1]);
        im8{1,3}(:,:,2) = imadjust(double(im8{1,3}(:,:,2)),[0 0.88],[0 1]);
        im8{1,4}(:,:,2) = imadjust(double(im8{1,4}(:,:,2)),[0 0.86],[0 1]);
        im8{1,5}(:,:,2) = imadjust(double(im8{1,5}(:,:,2)),[0 0.86],[0 1]);
        im8{1,6}(:,:,2) = imadjust(double(im8{1,6}(:,:,2)),[0 0.80],[0 1]);
        im8{1,7}(:,:,2) = imadjust(double(im8{1,7}(:,:,2)),[0 0.835],[0 1]); 

        im8{1,1}(:,:,3) = imadjust(double(im8{1,1}(:,:,3)),[0 0.75],[0 1]);
        im8{1,2}(:,:,3) = imadjust(double(im8{1,2}(:,:,3)),[0 0.85],[0 1]);
        im8{1,3}(:,:,3) = imadjust(double(im8{1,3}(:,:,3)),[0 0.83],[0 1]);
        im8{1,4}(:,:,3) = imadjust(double(im8{1,4}(:,:,3)),[0 0.815],[0 1]);
        im8{1,5}(:,:,3) = imadjust(double(im8{1,5}(:,:,3)),[0 0.79],[0 1]);
        im8{1,6}(:,:,3) = imadjust(double(im8{1,6}(:,:,3)),[0 0.74],[0 1]);
        im8{1,7}(:,:,3) = imadjust(double(im8{1,7}(:,:,3)),[0 0.765],[0 1]);       
    elseif (Expnr == 24 && t>8571)
        im8{1,1}(:,:,1) = imadjust(double(im8{1,1}(:,:,1)),[0 0.81],[0 1]);
        im8{1,2}(:,:,1) = imadjust(double(im8{1,2}(:,:,1)),[0 0.93],[0 1]);
        im8{1,3}(:,:,1) = imadjust(double(im8{1,3}(:,:,1)),[0 0.91],[0 1]);
        im8{1,4}(:,:,1) = imadjust(double(im8{1,4}(:,:,1)),[0 0.93],[0 1]);
        im8{1,5}(:,:,1) = imadjust(double(im8{1,5}(:,:,1)),[0 0.93],[0 1]);
        im8{1,6}(:,:,1) = imadjust(double(im8{1,6}(:,:,1)),[0 0.89],[0 1]);
        im8{1,7}(:,:,1) = imadjust(double(im8{1,7}(:,:,1)),[0 0.92],[0 1]);   

        im8{1,1}(:,:,2) = imadjust(double(im8{1,1}(:,:,2)),[0 0.80],[0 1]);
        im8{1,2}(:,:,2) = imadjust(double(im8{1,2}(:,:,2)),[0 0.93],[0 1]);
        im8{1,3}(:,:,2) = imadjust(double(im8{1,3}(:,:,2)),[0 0.92],[0 1]);
        im8{1,4}(:,:,2) = imadjust(double(im8{1,4}(:,:,2)),[0 0.95],[0 1]);
        im8{1,5}(:,:,2) = imadjust(double(im8{1,5}(:,:,2)),[0 0.95],[0 1]);
        im8{1,6}(:,:,2) = imadjust(double(im8{1,6}(:,:,2)),[0 0.90],[0 1]);
        im8{1,7}(:,:,2) = imadjust(double(im8{1,7}(:,:,2)),[0 0.86],[0 1]); 

        im8{1,1}(:,:,3) = imadjust(double(im8{1,1}(:,:,3)),[0 0.80],[0 1]);
        im8{1,2}(:,:,3) = imadjust(double(im8{1,2}(:,:,3)),[0 0.98],[0 1]);
        im8{1,3}(:,:,3) = imadjust(double(im8{1,3}(:,:,3)),[0 0.98],[0 1]);
        im8{1,4}(:,:,3) = imadjust(double(im8{1,4}(:,:,3)),[0 1],[0 1]);
        im8{1,5}(:,:,3) = imadjust(double(im8{1,5}(:,:,3)),[0 1],[0 1]);
        im8{1,6}(:,:,3) = imadjust(double(im8{1,6}(:,:,3)),[0 0.95],[0 1]);
        im8{1,7}(:,:,3) = imadjust(double(im8{1,7}(:,:,3)),[0 0.87],[0 1]);       
    else
        im8{1,1}(:,:,1) = imadjust(double(im8{1,1}(:,:,1)),[0 0.9],[0 1]);
        im8{1,2}(:,:,1) = imadjust(double(im8{1,2}(:,:,1)),[0 1],[0 1]);
        im8{1,3}(:,:,1) = imadjust(double(im8{1,3}(:,:,1)),[0 0.99],[0 1]);
        im8{1,4}(:,:,1) = imadjust(double(im8{1,4}(:,:,1)),[0 0.95],[0 1]);
        im8{1,5}(:,:,1) = imadjust(double(im8{1,5}(:,:,1)),[0 0.93],[0 1]);
        im8{1,6}(:,:,1) = imadjust(double(im8{1,6}(:,:,1)),[0 0.89],[0 1]);
        im8{1,7}(:,:,1) = imadjust(double(im8{1,7}(:,:,1)),[0 0.87],[0 1]);   

        im8{1,1}(:,:,2) = imadjust(double(im8{1,1}(:,:,2)),[0 0.86],[0 1]);
        im8{1,2}(:,:,2) = imadjust(double(im8{1,2}(:,:,2)),[0 0.98],[0 1]);
        im8{1,3}(:,:,2) = imadjust(double(im8{1,3}(:,:,2)),[0 0.92],[0 1]);
        im8{1,4}(:,:,2) = imadjust(double(im8{1,4}(:,:,2)),[0 0.95],[0 1]);
        im8{1,5}(:,:,2) = imadjust(double(im8{1,5}(:,:,2)),[0 0.97],[0 1]);
        im8{1,6}(:,:,2) = imadjust(double(im8{1,6}(:,:,2)),[0 0.92],[0 1]);
        im8{1,7}(:,:,2) = imadjust(double(im8{1,7}(:,:,2)),[0 0.96],[0 1]); 

        im8{1,1}(:,:,3) = imadjust(double(im8{1,1}(:,:,3)),[0 0.84],[0 1]);
        im8{1,2}(:,:,3) = imadjust(double(im8{1,2}(:,:,3)),[0 0.97],[0 1]);
        im8{1,3}(:,:,3) = imadjust(double(im8{1,3}(:,:,3)),[0 0.92],[0 1]);
        im8{1,4}(:,:,3) = imadjust(double(im8{1,4}(:,:,3)),[0 0.93],[0 1]);
        im8{1,5}(:,:,3) = imadjust(double(im8{1,5}(:,:,3)),[0 0.92],[0 1]);
        im8{1,6}(:,:,3) = imadjust(double(im8{1,6}(:,:,3)),[0 0.87],[0 1]);
        im8{1,7}(:,:,3) = imadjust(double(im8{1,7}(:,:,3)),[0 0.92],[0 1]);       
    end
    
    %% stitching
    if Expnr>20
        foto = [squeeze(im8{1,1}(745:3650,700:3470,:)) ...
        squeeze(im8{1,2}(836:3741,950:3640,:)) ...
        squeeze(im8{1,3}(888:3793,916:3780,:)) ...
        squeeze(im8{1,4}(785:3690,850:3820,:)) ...
        squeeze(im8{1,5}(835:3740,898:3450,:)) ...
        squeeze(im8{1,6}(960:3865,962:3680,:)) ...
        squeeze(im8{1,7}(750:3655,727:3770,:)) ...
        ];
    elseif Expnr>18
        foto = [squeeze(im8{1,1}(785:3690,700:3540,:)) ...
        squeeze(im8{1,2}(856:3761,850:3520,:)) ...
        squeeze(im8{1,3}(801:3706,946:3740,:)) ...
        squeeze(im8{1,4}(785:3690,850:3820,:)) ...
        squeeze(im8{1,5}(835:3740,898:3450,:)) ...
        squeeze(im8{1,6}(960:3865,962:3680,:)) ...
        squeeze(im8{1,7}(750:3655,727:3770,:)) ...
        ];
    elseif Expnr>10
        foto = [squeeze(im8{1,1}(745:3650,700:3470,:)) ...
        squeeze(im8{1,2}(776:3681,950:3640,:)) ...
        squeeze(im8{1,3}(768:3673,916:3780,:)) ...
        squeeze(im8{1,4}(785:3690,850:3820,:)) ...
        squeeze(im8{1,5}(835:3740,898:3450,:)) ...
        squeeze(im8{1,6}(960:3865,962:3680,:)) ...
        squeeze(im8{1,7}(750:3655,727:3770,:)) ...
        ];
    else    
        foto = [squeeze(im8{1,1}(745:3650,700:3470,:)) ...
        squeeze(im8{1,2}(776:3681,950:3640,:)) ...
        squeeze(im8{1,3}(821:3726,876:3740,:)) ...
        squeeze(im8{1,4}(785:3690,850:3820,:)) ...
        squeeze(im8{1,5}(835:3740,898:3450,:)) ...
        squeeze(im8{1,6}(960:3865,962:3680,:)) ...
        squeeze(im8{1,7}(750:3655,727:3770,:)) ...
        ];
    end  
    
    if ~(Expnr == 24 && t>8571)
        foto = imadjust(foto,[0.02 0.5],[0 1]);
    else
    end
%     figure(2);imshow(foto);
%     figure(5);imshow(foto(:,:,1));
%     figure(6);imshow(foto(:,:,2));
%     figure(7);imshow(foto(:,:,3));
    
%     %test = imfuse(im7{1,1}(745:3650,3602:3797,:),im7{1,2}(776:3681,1081:1276,:));
%     avg12 = imfuse(im8{1,1}(745:3650,3602:3797,:),im8{1,2}(776:3681,1081:1276,:),'blend');
%     %test = imfuse(im7{1,2}(776:3681,3620:3830,:),im7{1,3}(825:3730,863:1073,:));
%     avg23 = imfuse(im8{1,2}(776:3681,3620:3830,:),im8{1,3}(825:3730,863:1073,:),'blend');
%     %test = imfuse(im7{1,3}(825:3730,3710:3894,:),im7{1,4}(781:3686,829:1013,:));
%     avg34 = imfuse(im8{1,3}(825:3730,3710:3894,:),im8{1,4}(781:3686,829:1013,:),'blend');
%     %test = imfuse(im7{1,4}(781:3686,3720:3880,:),im7{1,5}(836:3741,810:970,:));
%     avg45 = imfuse(im8{1,4}(781:3686,3712:3881,:),im8{1,5}(835:3740,791:960,:),'blend');
%     %test = imfuse(im7{1,5}(835:3740,3310:3590,:),im7{1,6}(960:3865,821:1101,:));
%     avg56 = imfuse(im8{1,5}(835:3740,3310:3590,:),im8{1,6}(960:3865,821:1101,:),'blend');
%     %test = imfuse(im7{1,6}(960:3865,3670:3970,:),im7{1,7}(750:3655,737:1037,:));
%     avg67 = imfuse(im8{1,6}(960:3865,3670:3970,:),im8{1,7}(750:3655,737:1037,:),'blend');
%     
%         foto = [squeeze(im8{1}(745:3650,700:3601,:)) ...
%                 squeeze(double(avg12(:,:,:))/255) ...
%                 squeeze(im8{2}(776:3681,1277:3619,:)) ...
%                 squeeze(double(avg23(:,:,:))/255) ...
%                 squeeze(im8{3}(821:3726,1074:3709,:)) ...
%                 squeeze(double(avg34(:,:,:))/255) ...
%                 squeeze(im8{4}(785:3690,1014:3711,:)) ...
%                 squeeze(double(avg45(:,:,:))/255) ...
%                 squeeze(im8{5}(835:3740,961:3309,:)) ...
%                 squeeze(double(avg56(:,:,:))/255) ...
%                 squeeze(im8{6}(960:3865,1102:3669,:)) ...
%                 squeeze(double(avg67(:,:,:))/255) ...
%                 squeeze(im8{7}(750:3655,1038:3770,:)) ...
%                 ];    
    clear im7 im9
    
    %% RGB plot & save
    if rgb_img == 1
        imwrite(foto, ['STITCH\Exp' num2str(Expnr,'%03i') '_RGB_' num2str(t,'%05i') '.tif']);
        %         foto1 = flip(rot90(foto));
        %         R = georefcells([0 20],[0 3],size(foto1));
        %         geotiffwrite(['STITCH\Exp' num2str(Expnr,'%03i') '_RGB_' num2str(t,'%04i') '.tif'],foto1,R);
    end
    
    if rgb_fig == 1 || movieon_rgb == 1
        fig1 = figure('units','normalized','outerposition',[0 0.2 1 0.41],'visible','on');
        iptsetpref('ImshowAxesVisible','on');
        h=axes('Units','pixels','Position',[30 30 2000*0.93 300*0.93]);
        I=imshow(foto,'Parent',h,'XData',[0 20000],'YData', [0 3000],'Border','tight');
        set(h,'XTick',0:1000:20000,'XTickLabel',0:1:20,'XMinorTick','on', ...
            'YTick',0:1000:3000,'YTickLabel',{'3m','2','1','0'},'box','on', ...
            'TickDir','in','YMinorTick','on','YTickLabelRotation',270);
        if ~isnan(icam(t,1))
            title(['Exp', num2str(Expnr,'%03i'), '   Date = ', datestr(datetime(x.datetime(icam(t,1)),'ConvertFrom','datenum')), '   Cycle = ', num2str(t,'%05i')]);
        else
            title(['Exp', num2str(Expnr,'%03i'), '   Cycle = ', num2str(t,'%05i')]);
        end
        
        countmov_rgb = countmov_rgb + 1;
        if (movieon_rgb == 1 && countmov_rgb == int_mov/int_img)
            writeVideo(vid,getframe(gcf));
            countmov_rgb = 0;
        end
        
        countfig_rgb = countfig_rgb + 1;
        if (rgb_fig == 1 && countfig_rgb == int_fig/int_img)
            print(fig1,['FIGS\Exp' num2str(Expnr,'%03i') '_RGB_text_' num2str(t,'%05i')],'-dpng','-r200');
            countfig_rgb = 0;
        end
        
        close all
        close all hidden
        close all force
    end
    
    %% LAB plot & save
    if lab_fig == 1 || movieon_lab == 1 || lab_img == 1
        foto2 = rgb2lab(foto,'WhitePoint','e');
        %figure();imshow(foto2(:,:,1),[0 100]);
        %figure();imshow(foto2(:,:,2),[-127 127]);
        %figure();imshow(foto2(:,:,3),[-127 127]);
        clear foto
        foto2(:,:,3) = (foto2(:,:,3)+127)/255;
        %figure();imshow(foto2(:,:,1));
        %figure();imshow(foto2(:,:,2));
        %figure();imshow(foto2(:,:,3));
        foto5 = wiener2(foto2(:,:,3),[8 8]); % for save original stretch
        %figure();imshow(foto5);
        foto3(:,:,3) = imadjust(foto2(:,:,3),[0.36 0.49],[0 1]);
        %figure();imshow(foto3(:,:,3));
        clear foto2   
        foto4 = wiener2(foto3(:,:,3),[8 8]); % for save visible figure and movie
        %figure();imshow(foto4);
        clear foto3
    end
    
    if lab_img == 1
        imwrite(foto5, ['STITCH\Exp' num2str(Expnr,'%03i') '_LAB_' num2str(t,'%05i') '.tif']);
        %         foto6 = flip(rot90(foto5));
        %         R = georefcells([-0.02 20.02],[0.06 2.95],size(foto6));
        %         geotiffwrite(['STITCH\Exp' num2str(Expnr,'%03i') '_LAB_' num2str(t,'%04i') '.tif'],foto6,R);
    end
    
    if (lab_fig == 1 || movieon_lab == 1)
        fig2 = figure('visible','off');
        iptsetpref('ImshowAxesVisible','on');
        imshow(foto4,'XData',[0 20000],'YData', [0 3000],'Border','tight');
        colormap(gca,mymap);caxis([0 1]);
        text(size(foto4,2)/16,size(foto4,1)/16,['Exp' num2str(Expnr,'%03i') '   Date = ' num2str(x.date(icam(t,4))) '   Time = ' num2str(round((x.time(icam(t,4)))/100000)) '   Cycle = ' num2str(x.IMnr(icam(t,4)),'%05i')]);
        %mymap_blue=colormap(gca)
        %save('R:\Metronome\experiments\matlab_lisanne\mymap_blue','mymap_blue');
        
        countmov_lab = countmov_lab + 1;
        if (movieon_lab == 1 && countmov_lab == int_mov/int_img)
            writeVideo(vid2,getframe);
            countmov_lab = 0;
        end
        
        countfig_lab = countfig_lab + 1;
        if (lab_fig ==1 && countfig_lab == int_fig/int_img)
            print(fig2,['FIGS\Exp' num2str(Expnr,'%03i') '_LAB_text_' num2str(t,'%05i')],'-dpng','-r200');
            countfig_lab = 0;
        end

        clear foto4
        close all
        close all hidden
        close all force
    end
    
    t
end

%% MOVIE
if movieon_rgb == 1
    close(gcf), close(vid);
end
if movieon_lab == 1
    close(gcf), close(vid2);
end

toc
