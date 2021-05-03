function [latrng,lonrng] = tilegeorange(x_tile_rng,y_tile_rng,z)
%TILEGEORANGE calculates the lat and lon limits of a rectangular
%             array of geographic map tiles.
%
%   This function provides the inclusive lat and lon range of the map tiles
%   ranges input as x_tile_rng and y_tile_rng for a zoom level of z. Note
%   that even if a single tile number is requested, the outputs will be 
%   [1x2] because the lat/lon limits provide the outer bounds of the tiles
%   numbers provided. 
%
% Equations based on:
% https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
%
%
%
%Inputs: 
% x_tile_rng    is a vector [1x2] containting the tile number range in the 
%               x-direction (longitude) [x_left , x_right]
% y_tile_rng    is a vector [1x2] containting the tile number range in the 
%               x-direction (latitude) [y_bot , y_top]
% 
% z             is the zoom level of the map tiles [1x1] in the range 1-16
%
%Outputs:
% latrng        is vector [1x2] of [latbot lattop] with values in the 
%               range -90 to 90 deg
% lonrng        is vector [1x2] of [latbot lattop] with values in the 
%               range -180 to 180
%
%
%
% Author: Michal Shafer
% Date: 2019-05-20

% Equations based on:
% https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames

%%

if length(x_tile_rng)==1
        x_tile_rng(2) = x_tile_rng;
end
if length(y_tile_rng)==1
        y_tile_rng(2) = y_tile_rng;
end

n = 2 ^ z;
lon_deg_left = x_tile_rng(1) / n * 360.0 - 180.0;
lat_rad = atan(sinh(pi * (1 - 2 * y_tile_rng(2) / n)));
lat_deg_top = lat_rad * 180.0 /pi;

lon_deg_right = (x_tile_rng(2)+1) / n * 360.0 - 180.0;
lat_rad = atan(sinh(pi * (1 - 2 * (y_tile_rng(1)+1) / n)));
lat_deg_bot = lat_rad * 180.0 /pi;

latrng = [lat_deg_bot lat_deg_top];
lonrng = [lon_deg_left lon_deg_right];

end

