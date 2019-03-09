function [latlon] = xy2latlon(latlon_0,xy)
%function [lat,lon] = xy2latlon(latlon_0,xy)
%XY2LATLON Converts a series of x-y positions to lats and lons based on the
%given zero position in lat lon

%INPUTS:
% latlon_0  1x2    the [latitude longitude] of the origin (0,0)
%   xy      nx2    a matrix of nx2 of x (col 1) and y (column 2) of the position to be
%                  converted. +X is +East, +Y is +North. 

%Other variables
%	lat_0 is the latitude in degrees of the origin (0,0)
%	lon_0 is the longitude in degrees of the origin (0,0)
%   x is a list of x positions (+East) to be converted
%   y is a list of y positions (+North) to be converted
%   lon is the list of longitudes corresponding to the x y points
%   lat is the list of latitudes corresponding to the x y points

%See 2018-08-31 in research log notebook

if size(xy,2)~=2
    error('The xy matrix must be size n x 2, with x in column 1 and y in column 2')
end

if numel(latlon_0)~=2
     error('The origin lat-lon position must have only two elements')
end

num_pos = size(xy,1);

x = xy(:,1);
y = xy(:,2);

lat_0 = latlon_0(1);
lon_0 = latlon_0(2);


d = sqrt(x.^2+y.^2);
theta = atan2(x,y);%with +N, opposite is x and adjacent is y

r_earth = 6371000; %Earth radius in m
%lats2 = lat_0+180/pi*d.*cos(theta)/r_earth
%lons2 = lon_0+180/pi*acos((cos(-d/r_earth)-sind(lat_0).*sind(lats2))./(cosd(lat_0).*cosd(lats2)))

d_lat = 180/pi*d.*cos(theta)/r_earth;
lat = lat_0+d_lat;

the_argument = (cos(-d/r_earth)-sind(lat_0).*sind(lat))./(cosd(lat_0).*cosd(lat));
%Occasionally, the acos of this argument caused complex results because of
%values slightly larger than 1, so we round here to eliminate that
%possiblitity. 
the_argument_round = round(the_argument,15);
d_lon_abs = 180/pi*acos(the_argument_round);
%great circle distance lists the delta sigma equation as 'absolute value'
%I have to do the sign of y business to get the proper directions
d_lon = sign(x).*d_lon_abs;
lon = lon_0+d_lon;

latlon = [lat,lon];

end

