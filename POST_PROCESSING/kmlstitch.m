function [text_combined] = kmlstitch(flnm_list)
%KMLSTITCH reads a series of .kml files and merges their contents into a
%char array that can then be written as part of a larger kml file.
%   This function receives an cell array of filenames of KML files that the
%   user want to combine. The function reads in the text
%   of each one and then removes the contents after the document name and 
%   before the end of the document. 

%INPUTS:
%flnm_list      nx1 cell array of .kml files to be merged (including the
%               .kml as part of the filename.)


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

%Create output string. 
text_combined = text_hold;

end

