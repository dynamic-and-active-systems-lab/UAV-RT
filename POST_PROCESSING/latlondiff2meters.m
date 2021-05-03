function [dx,dy] = latlondiff2meters(latlon_1,latlon_2)
%latlondiff2meters Finds the distance in meters between two lat/lon
%positions. Inputs must have latitudes in the first column and longitudes
%in the second. If only one latlon entry is entered for latlon_1, all
%differences are simply made to this point, and visa-versa. 

%You can use atan2(dy,dx) to get the angle(rads) to between the two
%using a standard mathematical CS for the angle, with 0 on the East vector
%and positive angles counterclockwise.

%You can use atan2(dx,dy) to get the a standard compass bearing
%angle from North (in radians)

%INPUTS:
    %IMPORTANT: ENTRIES MUST BE ON DIFFERENT ROWS
% latlon_1  nx2    the [latitude longitude] of the origin (0,0)
% latlon_2  nx2    the [latitude longitude] of the point to measure to (0,0)

%OUTPUT:
%   dxy     nx2   a matrix of nx2 of x (col 1) and y (column 2) distance 
%                  from the first point to the second.
%                  +X is +East, +Y is +North. 


lat_1 = latlon_1(:,1);
lon_1 = latlon_1(:,2);
lat_2 = latlon_2(:,1);
lon_2 = latlon_2(:,2);

lats = [lat_1;lat_2];
lons = [lon_1;lon_2];

%Adapated from example https://www.mathworks.com/help/map/ref/mfwdtran.html
dczone = utmzone(mean(lats,'omitnan'),mean(lons,'omitnan'));
utmstruct = defaultm('utm'); 
utmstruct.zone = dczone;  
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);
[x1,y1] = mfwdtran(utmstruct,latlon_1(:,1),latlon_1(:,2));
[x2,y2] = mfwdtran(utmstruct,latlon_2(:,1),latlon_2(:,2));

dx = x2-x1;
dy = y2-y1;
end

