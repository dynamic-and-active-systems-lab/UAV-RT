function [flnmlist_mat] = downloadmaps(latlim,lonlim,zoomval,maptype,savedir)
%DOWNLOADMAPS downloads map tiles in a lat-lon range from the USGS WMS.
%
%   This function downloads map tiles in the range specified by the latlim
%   and lonlim inputs from the USGS Web Map Server at the zoom value
%   specified by 'zoom.' The maptype input string specifies which map
%   should be downloaded and the savedir string specifies where the maps
%   should be saved. The output is a matrix in the size of the tiles
%   requested that contains the filenames of the downloaded map tiles. Map
%   tiles are downloaded as 256x256 pixel .png files. 
%   
%   The files will be saved within the folder specified by savedir. It is
%   the caller's responsibility to generate these folders and setup a file
%   structure for different zoom levels, if desired. 
%
%   Before downloading, the function checks how many maps actually need to
%   be downloaded and checks with the user to make sure it is okay to
%   process. This just incase a large lat/lon limit and high zoom level
%   don't accidentally start a very large download request. 
%
%   The WMS servers often don't respond, and the download times out. The
%   timeout time is set to 3 second. When a download times out, the
%   funciton will throw a warning in the command line but will continue to
%   attempt the download five times. It usually only takes one more
%   download request, but if the server is down, 
%
%
%   The map names are give map_X***_Y***.png where *** indicates the tile
%   number. Note that the number of digits of *** will depend on the zoom
%   level and the lat/lon of that particular map. The digits are not fixed.


%Inputs:
%   latlim      Map latitude limits  [lat_S lat_N] (1x2)
%   lonlim      Map longitude limits [lon_W lon_E] (1x2)
%   zoom        Map tile zoom level 1-6 (1x1)
%   maptype     String of the map type to be downloaded. Valide inputs are 
%               'USGSHydroCached'
%               'USGSImageryOnly'
%               'USGSImageryTopo'
%               'USGSShadedReliefOnly'
%               'USGSTNMBlank'
%               'USGSTopo'
%               These can be found here: 
%               https://basemap.nationalmap.gov/arcgis/rest/services
%               The values must be exact because they are used to generate
%               urls that are used to download the .png files. 
%   savedir     String of the path (absolute or relative) to the directory
%               where the map cache should be saved. 


%Outputs:
%   flnmlist_mat    This is a cell array of the output file names in an nxm
%                   array, where the n rows correspond to latitudinal
%                   variation of the tiles and the m rows correspond the
%                   longitudinal variation of the tiles. 


%   Map tiles numbering uses the following coordinates in relation to the
%   lat,lon. More information can be found here:
%   https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames#Pseudo-code
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
%Author: Michael Shafer
%Date:   2019-05-22

%%
skip=false;%Flag used to skip the rest of the function if something goes wrong

%% CHECK OUT THE TOTAL NUMBER OF MAPS NEEDED. 
%First calculate the number of maps that are going to need to be
%downloaded, and check with the user to make sure this is okay.
num_of_maps = 0;
z = zoomval;
[x_tile_list,y_tile_list] = tilelookup(latlim,lonlim,z);
x = [x_tile_list(1), x_tile_list(end)];
y = [y_tile_list(1), y_tile_list(end)];
num_of_maps = (abs(diff(x))+1)*(abs(diff(y))+1)+num_of_maps;
approx_map_size = num_of_maps*30/1000;%Looks like all tiles are less than 30 kB. Report in Mb.
question = ['Are you sure you want to process with zoom level ',num2str(zoomval),' download? This will download ',num2str(num_of_maps),' and will require approximately ',num2str(approx_map_size),' Mb. Maps already downloaded will be skipped'];
answer = questdlg(question,'Yes','No');
if strcmp(answer,'Yes')
else %If either 'no' or 'cancel' are returned
    skip = true;
end

%% BEGIN MAP DOWNLOADS
if ~skip
    websaveoptions = weboptions('Timeout',3);%set so that downloads timeout after 3 s. WILL RETRY AFTER THAT.
    count = 0;
    f_waitbar = waitbar(count/num_of_maps,[maptype,': Dowloading zoom level XX maps.']);
    z = zoomval;
    [x_tile_list,y_tile_list] = tilelookup(latlim,lonlim,z);%Look up the required tile numbers
    %Server url definition
    baseurl = 'https://basemap.nationalmap.gov/arcgis/rest/services/${map}/MapServer/tile/${z}/${y}/${x}';
    
    try
        [~, msg, ~] = mkdir('mapcache');
    catch % In case makedirectory fails, shoot error as a warning
        warning(['Create of mapcache folder failed because: ',msg]);
        skip = true;
    end
    try
        [~, msg, ~] = mkdir(['mapcache/',num2str(z)]);
    catch % In case makedirectory fails, shoot error as a warning
        warning(['Create of map folders failed because: ',msg])
        skip = true;
    end
    timeout_tick = 0;
    if ~skip
        flnmlist = [];
        for j = 1:length(y_tile_list)
            for i = 1:length(x_tile_list)
                %Use the base URL and replace the map and xyz values to
                %generate the correct URL. 
                
                url_m = strrep(baseurl,'${map}',maptype);%Set the right map
                url_mz = strrep(url_m,'${z}',num2str(z));%Set the zoom level
                url_mzy = strrep(url_mz,'${y}',num2str(y_tile_list(j)));%Set the y tile
                url_mzyx = strrep(url_mzy,'${x}',num2str(x_tile_list(i)));%Set the x tile
                
                savename = ['map_x',num2str(x_tile_list(i)),'_y',num2str(y_tile_list(j)),'.png'];
                savepath = [savedir,'/',savename];
                if ~isfile(savepath)%Only download if the file doesn't already exist.
                    
                    try
                        websave(savepath,url_mzyx,websaveoptions);
                        timeout_tick = 10; %If websave doesn't fail this time, set to 10 so we know it was successful.
                    catch ME
                        warning(['Attempting websave again: ',ME.message])
                        timeout_tick = 1;
                        %Sometimes the server throws a
                        %502 error, but if you try
                        %again it works find
                        while strcmp(ME.identifier,'MATLAB:webservices:HTTP502StatusCodeError') && (timeout_tick<5)
                            try
                               websave(savepath,url_mzyx,websaveoptions);
                               timeout_tick = 10; %If websave doesn't fail this time, over set the timeout tick to get out of the loop.
                            catch ME
                                warning(['Attempting websave again: ',ME.message])
                                timeout_tick = timeout_tick+1;
                            end
                        end
                        if timeout_tick==5% if equal to 10, download was successful.
                            warning(['The map ',savestr,' was skippped because 5 download attemptes timed out.'])
                        end
                    end
                end
                if timeout_tick == 5
                    flnmlist_mat{j,i} = 'Not downloaded due to multiple timeouts.';
                else
                    flnmlist_mat{j,i} = savepath;
                end
                timeout_tick = 0;%Reset tick. 
                count = count+1;
                waitbar(count/num_of_maps,f_waitbar,[maptype,': Dowloading zoom level ',num2str(z),' maps.']);
            end
        end
        
    end
end
close(f_waitbar)
end

