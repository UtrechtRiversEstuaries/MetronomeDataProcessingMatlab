% This script creates a basic plot from processed Metronome DEM data in a netCDF file
% Written by Jan-Eike Rossius in September 2022

close all
clear

%% Input - to be filled in individually
% your personal working directory in which you use the standard folder structure:
pwd = 'C:\Users\7062052\OneDrive - Universiteit Utrecht\Open Science Metronome\MetronomeWorkingDirectory\';
% the current experiment (three digit number, in case of pilots use e.g. '052\Pilot1' as this string is used for directory paths)
exp = '052\Pilot4';
% define the name of the variable from the file that should be plotted
plotvar='z50';

%% Initialisation
%% Genaral

expalt=strrep(exp,'\','-'); % replacing \ in exp string to have a string that is nicer for naming files etc.

% open the file
fileid = netcdf.open([pwd '01metronome_experiments\Exp' exp '\processed_data\DEMs\laser_scanner\Exp' expalt 'DEMs.nc'],'NC_NOWRITE'); % open file for reading

% Get information about the contents of the file.
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(fileid);

% Identify location where DEMs are stored
if numdims==0 % If parent group of file doesn't contain any dimensions, search in groups
    childGrps = netcdf.inqGrps(fileid);
    for i=1:length(childGrps) % loop through groups
        groupName = netcdf.inqGrpName(childGrps(i));
        if strcmp(groupName,'DEMs') % find group called 'DEMs'
            DEMgroup=childGrps(i);
            [numDEMdims, numDEMvars, numDEMglobalatts, unlimDEMdimID] = netcdf.inq(DEMgroup); % get info about the DEM group
        end
    end
else % if no DEM group exists, the parent group acts as DEM group
    DEMgroup=fileid;
    [numDEMdims, numDEMvars, numDEMglobalatts, unlimDEMdimID] = netcdf.inq(DEMgroup);
end

%% Getting dimension data
for j=0:numDEMdims-1 % loop through dimensions to identify which dimension is which
    try
        axisAtt = netcdf.getAtt(DEMgroup,j,'axis'); % get the axis attribute
    catch exception
        axisAtt = 'none'; % put 'none' if axis attribute doesn't exist
    end
    if strcmp(axisAtt,'X') % if axis attribute is 'X', the current dimension is the x-dimension
        xdimid=j; % define the dimension id
        xdim=netcdf.getVar(DEMgroup,xdimid); % get the actual data contained
    elseif strcmp(axisAtt,'Y') % if axis attribute is 'Y', the current dimension is the y-dimension
        ydimid=j; % define the dimension id
        ydim=netcdf.getVar(DEMgroup,ydimid); % get the actual data contained
    elseif strcmp(axisAtt,'T') % if axis attribute is 'T', the current dimension is the time dimension
        tdimid=j; % define the dimension id
        tdim=netcdf.getVar(DEMgroup,tdimid); % get the actual data contained
    end
end
% check whether all dimensions exist
existtest = exist('xdim','var');
existtest = existtest + exist('ydim','var');
existtest = existtest + exist('tdim','var');
% throw errors if not all three dimensions are found
if existtest==0
    error('All dimesions do not contain axis attribute');
elseif existtest<3 && existtest>0
    error('Only one or two dimensions contain axis attributes. It should be three');
end

%% Get variable data
varids = netcdf.inqVarIDs(DEMgroup); % get the variable ids
for k=1:numDEMvars % loop through variables
    [varname,xtype,dimids,natts] = netcdf.inqVar(DEMgroup,varids(k)); % get info bout the variable
    if strcmp(plotvar,varname)
        assignin('base','plotvar',netcdf.getVar(DEMgroup,varids(k))); % if not, save under the name of the variable
        plotvarid=varids(k);
    end
end

%% Determining figure and plot size
% get information about the screen size
set(0,'units','pixels')
screenpixels = get(0,'screensize');
figheight=screenpixels(4)-150; % set height of the figure to 150 pixels less than the height of the screen
numplotsvert=ceil(length(tdim)/2); % number of rows of plots
% check whether there is an even or odd number of plots necessary
if numplotsvert==length(tdim)/2
    even=true;
else
    even=false;
end

padding=75; % space at the sides of the figure used for axis labels etc.
vertpixplot=(figheight-padding)/(numplotsvert+1*double(even)); % determine number of vertical pixels per subplot. padding is left for axis description. Space of one plot row is left for colourbar if there's an even number of plots necessary.

if netcdf.getAtt(DEMgroup,xdimid,'units')~=netcdf.getAtt(DEMgroup,ydimid,'units') % check if units of both spatial axes are equal
    error('x- and y-dimensions do not have the same units. Adjust code manually!') % if not code needs to be adjusted as this is ssumed in the fllowing section
else
    pixperunit=vertpixplot/(max(ydim)-min(ydim)); % caculate number of pixels available per unit length to determine width of plots and figure
end

horzpixplot=(max(xdim)-min(xdim))*pixperunit; % calculate number of horizontal pixels per plot
figwidth=padding+2*horzpixplot; % set figure width to twice the plot width plus padding for axis labels

if figwidth>screenpixels(3) % in case the figure gets too wide, base the figure size on the screen width instead
    figwidth=screenpixels(3)-50;
    horzpixplot=(figwidth-padding)/2;
    pixperunit=horzpixplot/(max(xdim)-min(xdim));
    vertpixplot=pixperunit*(max(ydim)-min(ydim));
    figheight=vertpixplot*(numplotsvert+1*double(even))+padding;
end

% Creating figure
figure1=figure('Name',netcdf.getAtt(DEMgroup,-1,'title'), 'Position',[10 50 figwidth figheight]);

% Plot size and padding relative to figure size
vertrelplot=vertpixplot/figheight;
horzrelplot=horzpixplot/figwidth;
vertrelpad=padding/figheight;
horzrelpad=padding/figwidth;

%% Actual plotting
valrange=prctile(plotvar,[2 12 17 20 25 35 40 60 75 85 98],'all')'; % percentile values for an adjusted colormap
demcolorpos=valrange-valrange(1); % setting lowest value to zero and adjust the others
demcolorpos=demcolorpos./demcolorpos(end); % norming all values on the range 0-1
demcolorcol=flipud([0.51372549	0.235294118	0.047058824; 0.749019608	0.560784314	0; 0.894117647	0.764705882	0.290196078; 0.57254902	0.815686275	0.31372549; 0	0.690196078	0.31372549; 0.215686275	0.337254902	0.137254902; 0.866666667	0.921568627	0.968627451; 0.596078431	0.752941176	0.894117647; 0	0.690196078	0.941176471; 0.121568627	0.305882353	0.470588235; 0	0.125490196	0.376470588]); % colors for the colormap, ranging from blue via white, green and yellow to brown
demcolormap=flip(customcolormap(demcolorpos,demcolorcol)); % make colormap from positions and colors previously defined
fontsize=10; % the fontsize used in the plots
fontheight=15; % pixels for one line of text
pixperlabel=25; % minimum number of pixels that should be available in one subplot for one axis label on the y-axis
numyticks=floor(vertpixplot/pixperlabel)+1; % calculating the number of y-ticks
for l=1:numyticks
    yticks(l)=round((min(ydim)+(l-0.5)*((max(ydim)-min(ydim))/numyticks))*100)/100; % determining the position of the y-ticks
end
numxticks=floor(horzpixplot/(2*pixperlabel))+1; % calculating the number of x-ticks (twice the minimum distance as for y-axis)
for l=1:numxticks
    xticks(l)=round((min(xdim)+(l-0.5)*((max(xdim)-min(xdim))/numxticks))*100)/100; % determining the postion of the x-ticks
end
for i=1:length(tdim) % looping through the time-dimension
    if i<=numplotsvert % define positions of plots in figure for the left column
        posx=horzrelpad*0.75;
        posy=1-vertrelpad*0.25-((i+double(even))*vertrelplot);
    else % define positions of plots in figure for the right column
        posx=horzrelpad*0.75+horzrelplot;
        posy=1-vertrelpad*0.25-((i-numplotsvert+1)*vertrelplot);
    end

    ax(i)=axes('position',[posx posy horzrelplot vertrelplot]); % creating subplot
    imagesc(xdim,ydim,plotvar(:,:,i)'); % plotting respective timestep
    colormap(ax(i), demcolormap); % set colormap
    text(min(xdim)+((max(xdim)-min(xdim))/100),max(ydim)-((fontheight/vertpixplot)*(max(ydim)-min(ydim))),num2str(tdim(i),'%04.0f'),'FontSize',fontsize); % label the value of the time step
    text(min(xdim)+((max(xdim)-min(xdim))/100),max(ydim)-2*((fontheight/vertpixplot)*(max(ydim)-min(ydim))),netcdf.getAtt(DEMgroup,tdimid,'units'),'FontSize',fontsize,'Interpreter','none'); % label the unit of the time step
    caxis(valrange([1 end])) % set the color axis
    set(gca,'YDir','normal'); % set y-axis direction to upward positive
    axis equal; % couple the axis scale, note that units should be equal
    set(gca, 'YTick', yticks,'FontSize',fontsize) % set y-ticks
    set(gca, 'XTick', xticks,'FontSize',fontsize) % set x-ticks

    if i>numplotsvert
        set(gca,'YTickLabels',[]); % no y-axis label for right column of plots
    else
        ylabel(['y [' netcdf.getAtt(DEMgroup,ydimid,'units') ']'],'FontSize',fontsize); % label the y-axis with corresponding unit
    end

    if i==numplotsvert || i==length(tdim) 
        xlabel(['x [' netcdf.getAtt(DEMgroup,xdimid,'units') ']'],'FontSize',fontsize); % label the x-axis with corresponding unit for bottom plots
    else
        set(gca,'XTickLabels',[]); % no x-axis label for plots not at the bottom
    end
end
c=colorbar; % create colrbar
c.Location='northoutside'; % colorbar not within the plots
if even % for two plot columns with equal number of plots
    cbarwidth=2*horzrelplot; % wide colorbar
    cbarxpos=horzrelpad*0.75; % edge at the same position of plots
else % for two plot columns with different number of plots
    cbarwidth=0.9*horzrelplot; % smaller colorbar
    cbarxpos=horzrelpad*0.75+1.1*horzrelplot; % in the position of the "missing" plot at the top right
end
c.Position=[cbarxpos 1-vertrelpad*0.25-vertrelplot*0.3 cbarwidth vertrelplot*0.3]; % putting together colorbar position
set(c, 'YAxisLocation','bottom') % axis label at the bottom of colorbar
c.Label.String=[netcdf.getAtt(DEMgroup,plotvarid,'long_name') ' [' netcdf.getAtt(DEMgroup,plotvarid,'units') ']']; % get label for the colorbar from variable
c.Label.FontSize=fontsize; % set label fontsize
c.FontSize=fontsize; % set ticks fontsize

%% Closing file
netcdf.close(fileid);

%% Functions
% ----------------------------------------------------------------------- %
% FUNCTION "customcolormap" defines a customized colobar given the        %
% positions and the colors that are going to generate the gradients.      %
%                                                                         %
%   Input parameters:                                                     %
%       - positions:    Vector of positions, from 0 to 1. Note that the   %
%                       first position must be 0, and the last one must   %
%                       be 1.                                             %
%       - colors:       Colors to place in each position. This parameter  %
%                       can be specified as a RGB matrix (n_colors x 3), or
%                       as a cell vector, containing HTML values.         %
%                       For instance: {'#ffffff','#ff0000','#000000'} is  %
%                       equivalent to [1 1 1; 1 0 0; 0 0 0].              %
%       - m:            (Optional) Number of points (recommended: m > 64).%
%                                                                         %
%   Output variables:                                                     %
%       - J:            Colormap in RGB values (dimensions [mx3]).        %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       J = customcolormap([0 0.5 1], {'#ffffff','#ff0000','#000000'});   %
%       colorbar; colormap(J);                                            %
% ----------------------------------------------------------------------- %
%   Versions:                                                             %
%       - v1.0.:    (19/11/2018) Original script.                         %
% ----------------------------------------------------------------------- %
%       - Author:   Víctor Martínez-Cagigal                               %
%       - Date:     19/11/2018                                            %
%       - Version:  1.0                                                   %
%       - E-mail:   vicmarcag (at) gmail (dot) com                        %
%                                                                         %
%       Biomedical Engineering Group (University of Valladolid), Spain    %
% ----------------------------------------------------------------------- %
function J = customcolormap(positions, colors, m)

    % Error detection and defaults
    if nargin < 3
       f = get(groot,'CurrentFigure');
       if isempty(f)
          m = size(get(groot,'DefaultFigureColormap'),1);
       else
          m = size(f.Colormap,1);
       end
    end
    if ~isnumeric(m), error('Parameter m must be numeric.'); end
    if nargin < 2, error('Not enough parameters.'); end
    
    if iscell(colors)
        colors = colors(:);
        n_colors = length(colors);
        for i = 1:n_colors
            temp = colors{i};
            if ~ischar(temp)
                error(['Colors must be specified in HEX format (e.g., #FFFFFF).' ...
                ' Type "help colorbar" for further information']);
            elseif ~strcmp(temp(1),'#')
                error(['Character # is missing if %s.' ...
                ' Type "help colorbar" for further information'], temp);
            elseif length(temp)~=7
                error(['Not a valid color format: %s (use this format: #FFFFFF).' ...
                ' Type "help colorbar" for further information'], temp);
            end
        end
    elseif ismatrix(colors)
        n_colors = size(colors);
        if length(n_colors) ~= 2
            error(['RGB colors must be a 2D matrix.' ...
                ' Type "help colorbar" for further information']);
        elseif n_colors(2) ~= 3
            error(['RGB colors matrix must have 3 columns.' ...
                ' Type "help colorbar" for further information']);
        elseif min(colors(:))<0 || max(colors(:))>255
            error(['RGB colors matrix values must range from 0 to 255.' ...
                ' Type "help colorbar" for further information']);
        end
    else
        error(['Colors must be a cell vector or a matrix of RGB values.' ...
            ' Type "help colorbar" for further information']);
    end
    if ~isvector(positions)
        error(['Positions must be specified as a vector.' ...
            ' Type "help colorbar" for further information']);
    elseif min(positions)<0 || max(positions)>1
        error(['Positions must range from 0 to 1 in an ascending order.' ...
            ' Type "help colorbar" for further information']);
    elseif length(positions) ~= length(unique(positions))
        error(['Check the positions vector, there are some duplicates.' ...
            ' Type "help colorbar" for further information']);
    else
        positions = sort(positions, 'ascend');
        if positions(1)~=0
            error(['The first positions must be 0.' ...
            ' Type "help colorbar" for further information']);
        elseif positions(length(positions))~=1
            error(['The last positions must be 1.' ...
            ' Type "help colorbar" for further information']);
        elseif length(positions) ~= n_colors
            error(['The number of positions does not match the number of colors.' ...
            ' Type "help colorbar" for further information']);
        end
    end

    % Convert HEX colors into RGB colors if required
    if iscell(colors)
        hex_colors = colors;
        colors = NaN(n_colors,3);
        for i = 1:n_colors
            colors(i,:) = hex2rgb(hex_colors{i});
        end
    end
    
    % Compute positions along the samples
    color_samples = round((m-1)*positions)+1;

    % Make the gradients among colors
    J = zeros(m,3);
    J(color_samples,:) = colors;
    diff_samples = diff(color_samples)-1;
    for d = 1:1:length(diff_samples)
        if diff_samples(d)~=0
            color1 = colors(d,:);
            color2 = colors(d+1,:);
            G = zeros(diff_samples(d),3);
            for idx_rgb = 1:3
                g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
                g([1, length(g)]) = [];
                G(:,idx_rgb) = g';
            end
            J(color_samples(d)+1:color_samples(d+1)-1,:) = G;
        end
    end
    J = flipud(J);
end

% Function that converts an HEX string (format: #FFFFFF) to the
% corresponding RGB color
function rgb = hex2rgb(hexString)
	if size(hexString,2) ~= 7
		error('Not a color! %s', hexString);
	else
		r = double(hex2dec(hexString(2:3)))/255;
		g = double(hex2dec(hexString(4:5)))/255;
		b = double(hex2dec(hexString(6:7)))/255;
		rgb = [r, g, b];
	end
end
