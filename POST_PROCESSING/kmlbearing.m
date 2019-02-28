function [] = kmlbearing(filename,lats,lons,alts,bear,drawdist,color,alpha_const)
%KMLBEARING creates a .kml file of bearning lines 
%   This function receives a series of latitudes and longitudes and draws a
%   series of lines eminating from those points a distance of drawdist at
%   a compass bearing (from true north) of bear. The color and transparancy
%   can be controlled as well. 

%INPUTS:
%filename       a char array of the filename to use for the kml
%lats           nx1 numeric vector of latitudes
%lons           nx1 numeric vector of longitudes
%alts           nx1 numeric vector of altitudes
%bear           nx1 numeric vector of compass bearings
%drawdist       nx1 numeric values of distance to draw line in meters
%color          1x3 RGB color vector
%alpha_conts    1x1 numeric values alpha values for transparency 0-1

%OUTPUTS:
%none except a KML file of filename.kml saved in the current directory. 


if ~isequal(numel(lats),numel(lons),numel(alts),numel(bear)) %Check to see that vectors are same length
    error('Number of latitude, longitude, altitude, and bearing elements must be equal')
else
    num_bears = length(lats);
end

%Make sure everything is the same size. 
if isrow(lats); lats = lats'; end
if isrow(lons); lons = lons'; end
if isrow(bear); bear = bear' ;end
if isrow(drawdist); drawdist = drawdist'; end
%if isrow(color); color = color'; end
if isrow(alpha_const); alpha_const = alpha_const' ;end


%Create a series of secondary points at distance from each original point
    r_earth = 6371000; %Earth radius in m
    lats2 = lats+180/pi.*drawdist.*cosd(bear)./r_earth;
    d_lons = 180/pi*acos((cos(-drawdist/r_earth)-sind(lats).*sind(lats2))./(cosd(lats).*cosd(lats2)));
    lons2 = lons+sign(sind(bear)).*d_lons;
    alts2 = alts.*ones(size(lats));

%Generate a series of temporary KML files that will later be stitched
%together. Each of these will be a single bearing line
kml_index = 1;
if strcmp(filename(1:end-3),'.kml')
    error('KML bearing filename must include the file extension .kml')
else
    temp_bear_path = filename(1:end-4);
end
mkdir(temp_bear_path)
for i = 1:num_bears
        flnm_lst_bear{kml_index} = [temp_bear_path,'/temp_bear_',num2str(kml_index),'.kml'];
        kmlwriteline(flnm_lst_bear{kml_index},[lats(kml_index) lats2(kml_index)],[lons(kml_index) lons2(kml_index)],[alts(kml_index) alts2(kml_index)],'AltitudeMode','relativeToGround','Color',color,'Name',['Bearing-',num2str(i)],'Alpha',alpha_const,'LineWidth', 1.5)
     kml_index = kml_index+1;
end


%Stitch all the temporary kml files together into a single character array
text_stitched = kmlstitch(flnm_lst_bear);
text_out = kmlwriteprep({text_stitched},'bearings');

%Write the stitched kml text to a compiled KML file
fid = fopen(filename,'wt');
fprintf(fid, text_out);
fclose(fid);

rmdir(temp_bear_path,'s')% remove temp directory and its contents
% %Delete the temporary kml files
% for i = 1:length(flnm_lst_bear)
%     delete(flnm_lst_bear{i})
% end


end

