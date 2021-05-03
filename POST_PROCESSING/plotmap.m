function [plottedmaps, wantedmaps] = plotmap(ax,latcent,loncent,latlim,lonlim,z,cachedir,plotselect)
%PLOTMAP plots downloaded map tiles in the range requested.
%
%   This function receives an axis (ax) as an input and plots a series of
%   map tiles on that axis. The function receives a bounds of lat and lon,
%   as well as the center point lat/lon, zoom level, the map cacher
%   directory on the local machine, and the plotselect. The fucntion first
%   determines which tile files are needed over the lat and lon limit
%   inputs. It then looks in the cache directory to see which maps are
%   available for plotting. The function then plots the maps that are
%   needed and avaialble in the cache directory. Each tile is plotted
%   individually. Tiles are 256 x 256 pixels. The x and y axes can either
%   be latitude and longitude if 'geo' is selected as the plotselect
%   variable, or in x and y distance in meters from the lat and lon center
%   inputs. If plotting in distance, the plotselect variable also include a
%   coordinate definitions with options of meters-NXEY or meters-NYEX
%   specify how North and East related to x and y. 


%Inputs:
%   ax          Matlab axis on which to plot
%   latcent     Map center Latitude  (1x1)
%   loncent     Map center Latitude  (1x1)
%   latlim      Map latitude limits  [lat_S lat_N] (1x2)
%   lonlim      Map longitude limits [lon_W lon_E] (1x2)
%   z           Map tile zoom level 1-16 (scalar 1x1)
%   cachedir    String of the path (absolute or relative) to the directory
%               where the map cache is located. File structure within this
%               directory must have folders named with the zoom level of
%               the tiles within that folder.
%   select      String selecting if plot should be in lat/lon (geo) or in
%               distance (meters-****). The **** term shoudl specify the
%               relationship betweeen the x and y axis (horizontal and
%               vetical, respectively) and the North and East directions.
%               For example, 'meters-NXEY' will produce a plot with the image
%               set up with North along the X axis and East along the Y.
%               Selecting 'meters-NYEX' will produce a plot with the image
%               set up with North along the Y axis and East along the X. In
%               either case, the image is plotted as a surface at z=0

%Outputs:
%   plottedmaps The number of maps actually plotted 
%   wantedmaps  The number of maps that would have been plotted if
%               available.

%Author: Michael Shafer
%Date:   2019-05-22

%%

hold(ax,'on');%Retain what is already on the plot
viewhold = ax.View;
%%Setup the slash character depending on the platform
if ismac
    slash_char = '/';
elseif ispc
    slash_char = '\';
end

[x_tile_list,y_tile_list] = tilelookup(latlim,lonlim,z);

wishlist = cell(length(y_tile_list),length(x_tile_list));
wishlist_path = cell(length(y_tile_list),length(x_tile_list));
%First make a list of all the 
for i = 1:length(x_tile_list)
    for j = 1:length(y_tile_list)
    savename = ['map_x',num2str(x_tile_list(i)),'_y',num2str(y_tile_list(j)),'.png'];
    wishlist{j,i} = savename;%A listing (in the same layout as should be plotted) of the tiles we'd like
    wishlist_path{j,i} = [cachedir,slash_char,num2str(z),slash_char,savename];%A listing (in the same layout as should be plotted) of the tiles we'd like
    end
end

%Find all the map files in the directory
allfilelist = dir([cachedir,slash_char,num2str(z)]);
tick = 1;
havelist = {};  %A list of all maps in the file
wishhavelist = {};%Set up as empty initially. Will fill in and explain later
for i = 1:length(allfilelist)%
    %First check to see if name is even long enough. Then check to see if
    %it has the name features of the cached maps
    if (length(allfilelist(i).name)>5) &&strcmp(allfilelist(i).name(1:5),'map_x') && strcmp(allfilelist(i).name(end-3:end),'.png') 
         havelist{tick} = allfilelist(i).name;
         havelist_path{tick} = [allfilelist(i).folder,allfilelist(i).name];
         tick=tick+1;
    end
end

if ~isempty(havelist)
%Now pull all the x and y tile numbers out of the file names
x_nums = zeros(1,length(havelist));
y_nums = zeros(1,length(havelist));
for i = 1:length(havelist)
    u_score_2_pos = find(havelist{i}=='_',1,'last');
    x_nums(i) = str2double(havelist{i}(6:u_score_2_pos-1));
    y_nums(i) = str2double(havelist{i}(u_score_2_pos+2:end-4));
end


%Here we plot each map file that we have downloaded that is also in the
%wishlist
wishhavelist = cell(numel(havelist),1);%A list of all maps in the file that we want to plot. Will be empty elements if the file is not in our wishlist
wishhavelist_path = cell(numel(havelist),1);%
for i = 1:length(havelist)
    found_ind = find(strcmp(wishlist,havelist{i}),1);
    if ~isempty(found_ind)
        wishhavelist{i} = wishlist{found_ind};%There will be empty matrices in the elemenets of wishhavelist that aren't in our wishlist
        wishhavelist_path{i} = wishlist_path{found_ind};
        
        [latrng,lonrng] = tilegeorange(x_nums(i),y_nums(i),z);
        
        dlat = latrng(2)-latrng(1);
        dlon = lonrng(2)-lonrng(1);
        
       
        image_curr = flip(imread(wishhavelist_path{i}),1);
                            %   Matlab uses the following coordinates for
                            %   the image display, so we flip the columns 
                            %   of the output to align with lat (+up) and
                            %   lon(+right) directions
                            %   about 
                            %  _______________> +x
                            %  |
                            %  |
                            %  |     IMAGE MATRIX (OUT)
                            %  |
                            %  |
                            %  V
                            %  +y
        %Build the map tag that we will tag to this plot so we can
        %reference them later in other codes if we need to.
        maptag = ['MAP-_x',num2str(x_nums(i)),'_y',num2str(y_nums(i)),'_z',num2str(z)];
        
        if strcmp(plotselect,'geo')
            %imshow(image_curr,'XData',lonrng,'YData',latrng,'Parent',ax,'Tag',maptag);
            %%The line above worked on a mack but not on a PC. Had to
            %%change it to set the tag manualy.
            imag_h = imshow(image_curr,'XData',lonrng,'YData',latrng,'Parent',ax);
            imag_h.Tag = maptag;
            hold(ax,'on');
        else
            [arclen_x_right,bear_x_right] = distance(latcent,loncent,latcent,lonrng(2));
            [arclen_x_left,bear_x_left]   = distance(latcent,loncent,latcent,lonrng(1));
            [arclen_y_top,bear_y_top]     = distance(latcent,loncent,latrng(2),loncent);
            [arclen_y_bot,bear_y_bot]     = distance(latcent,loncent,latrng(1),loncent);
            
            %Useful debugging 
            %deg2km(arclen_x_right-arclen_x_left)
            %deg2km(arclen_y_bot-arclen_y_top)
            
            km_x_right = deg2km(arclen_x_right)*sind(wrapTo180(round(bear_x_right)));
            km_x_left = deg2km(arclen_x_left)*sind(wrapTo180(round(bear_x_left)));
            km_y_top = deg2km(arclen_y_top)*cosd(wrapTo180(round(bear_y_top)));
            km_y_bot = deg2km(arclen_y_bot)*cosd(wrapTo180(round(bear_y_bot)));
            
            x_dist_data = 1000*linspace(km_x_left,km_x_right,size(image_curr,2));
            y_dist_data = 1000*linspace(km_y_bot,km_y_top,size(image_curr,1));

            %imshow(image_curr,'XData',x_dist_data,'YData',y_dist_data,'Parent',ax);
            [X_DIST_MAT, Y_DIST_MAT] = meshgrid(x_dist_data,y_dist_data);
            if strcmp(plotselect,'meters-NXEY')
                surface(ax,Y_DIST_MAT, X_DIST_MAT, zeros(size(X_DIST_MAT)), image_curr, 'facecolor', 'texturemap', 'edgecolor', 'none','Tag',maptag);
            elseif strcmp(plotselect,'meters-NYEX')
                surface(ax,X_DIST_MAT, Y_DIST_MAT, zeros(size(X_DIST_MAT)), image_curr, 'facecolor', 'texturemap', 'edgecolor', 'none','Tag',maptag);
            end
            hold(ax,'on');
        end
        %Display for debugging
        %disp('I plotted a map!')
    else
       %Only plot if we have a map
    end
end

if strcmp(plotselect,'geo')
        %set(ax,'YTick',y_lat_data);
        set(ax,'YDir','normal');
        %set(ax,'XTick',x_lon_data);
        set(ax,'XTickMode','auto')
        set(ax,'YTickMode','auto')
        set(ax,'Visible','on')
         set(ax,'XLim',lonlim)
         set(ax,'YLim',latlim)
        set(ax,'DataAspectRatio',[[1,dlat/dlon*256/256],1]); %[[1,dlat/dlon*size(out,2)/size(out,1)],1])
        hold(ax,'on');
%Set aspect ratio to scale by number of pixels
%and then by the delta lat/lon ratio. If one
%tile, the pixels would be 256x256 and would
%need a 1x1 aspect ratio, but because the x and
%y data are in lat and lon, the figure would be
%skewed because if you just used 'axis equal,'
%because lat and lon increments aren't the same
%scale (distance). You have to scale the
%splotting by the ratio of the delta lat and
%lon to maintain the 1x1 aspect ration of each
%pixel. 
        
elseif strcmp(plotselect,'meters')
        set(ax,'YDir','normal');
        set(ax,'XTickMode','auto')
        set(ax,'YTickMode','auto')
        set(ax,'Visible','on')
        %Don't need to change aspect ratio when plotting in meters, because
        %the x and y data has a unity aspect ratio. 
        hold(ax,'on');
end
ax.View = viewhold;
    
end
    plottedmaps = sum(~cellfun(@isempty,wishhavelist)); %Let the call know how many maps were actually plotted
    wantedmaps  = numel(wishlist); %Let the call know how many maps were actually plotted
end

