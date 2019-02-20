function [text_combined] = kmlfolder(flnm_list,foldername)
%KMLFOLDER reads a series of .kml files and merges their contents into a
%kml folder.
%   This function receives an cell array of filenames of KML files that the
%   user want to combine into a kml folder. The function reads in the text
%   of each one and then pulls the contents after the document name. All 
%   of these are then compiled into one char array and appended at the 
%   beginning and end with <Folder> <name> foldername</name> and </Folder>,
%   respectively. kmlwriteprep can then be used to complete the file and
%   make a valide kml text that can be written to a file. 

%INPUTS:
%flnm_list      nx1 cell array of .kml files to be merged (including the
%               .kml as part of the filename.)
%foldername     char array of the folder name the kmls will be grouped into

%OUTPUTS:
%a char array that contains all the text of the combined kml files that can
%then be writted to a .kml file. 


flnm_in = flnm_list;

%The number of files to read
num_fls = length(flnm_in);

%Find the place in the first file where the document termination string
%occurs
text_hold = '';
for i = 1:num_fls
    text_i_select = kmlextract(flnm_in{i});
    text_hold = [text_hold,text_i_select];%Add it to the compiled text string. 
end

NEW_KML_HEAD = sprintf(['\t<Folder>\r\n\t<name>',foldername,'</name>']);
NEW_KML_TAIL = sprintf(['\r\n\t</Folder>\r\n']);


%Merge the first file with all the others and stich on the document
%termination text to make it a valid kml. 
text_combined = [NEW_KML_HEAD,text_hold,NEW_KML_TAIL ];

end

