function [x_tile_list,y_tile_list] = tilelookup(latlim,lonlim,z)
%TILELOOKUP gets the x and y tile numbers of over a lat and lon range.
%
%   This function outputs the an x and y map tile list that covers the
%   lat and lon limits requested. A single lat and lon can be entered, 
%   in which case a single tile number will be output. If a range is
%   provided, a range of x and y tile numbers will be ouptut. If the
%   lat and/or lim range is less than the size of the tile, there will
%   only be a single tile number provided for x or y. 
%
% Equations based on:
% https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
%
%Inputs: 
% lonlim    is longitude limit [lonwest, loneast] in the range of -180 to
%           180. One value can be entered if only one lat/lon is needed. 
%           Can be [1x1] or [1x2]
% latlim    is latiude limit [latsouth, latnorth] in the range of -90 to
%           90. One value can be entered if only one lat/lon is needed. 
%           Can be [1x1] or [1x2]
% z         is zoom level (1-16) [1x1]
%
%
%Outputs: 
% x_tile_list  a 1xm list of x tile numbers that encompass the longitude
%              range requested. 
% y_tile_list  a 1xn list of y tile numbers that encompass the latitude
%              range requested. 
%
%   Map tiles numbering uses the following coordinates in relation to the
%   lat,lon. More information can be found here:
%   https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
%  
%  +90,-180_______________> +x       +90,+180
%       |
%       |
%       |     TILE NUMBER (x,y)
%       |
%       |
%       V
%       +y
%
%  -90,-180                          -90,+180
%
%
%
%
% Author: Michal Shafer
% Date: 2019-05-20

%%
if size(latlim)~=size(lonlim)
    error('Latlim and lonlim sizes must be equal')
end

max_x_tile = 2^z;%Maximum num of tiles in the x and  directions
max_y_tile = 2^z;
    %Get the first and last needed x tiles ([1x2])
    x = floor((lonlim+180)/360*2.^z);
    %Get the first and last needed y tiles ([1x2])
	y = floor((1-(log(tan(latlim.*pi/180)+1./(cos(latlim.*pi/180))))./pi).*2.^(z-1));
    
    if length(x)>1 
        %Generate a complete x tile list 
        x_tile_list = sort(linspace(x(1),x(2),abs(diff(x(:)))+1));
    else %If only one tile is needed the entire list is the output above
        x_tile_list = x;
    end
    if length(y)>1  
        y_tile_list = sort(linspace(y(1),y(2),abs(diff(y(:)))+1));
    else
        y_tile_list = y;
    end
    
    x_tile_list = wrapToN(x_tile_list,max_x_tile);
    y_tile_list = wrapToN(y_tile_list,max_y_tile);
    
    
    function out = wrapToN(in,n)
        %This function wraps a list on the range of 0-N. It interps to 360,
        %then uses matlab wrapTo360 function, then uninterps and sorts. 
        inon360 = interp1([0 n],[0 360],in,'linear','extrap');
        outon360 = wrapTo360(inon360);
        out = sort(interp1([0 360],[0 n],outon360));
        
    end
end

