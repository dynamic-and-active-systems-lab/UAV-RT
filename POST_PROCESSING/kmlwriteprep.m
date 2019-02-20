function [text_combined] = kmlwriteprep(charinputs,docname)
%KMLWRITEPRE combines a series of text inputs and adds the necessary 
%characters to the beginning and end so that it can be written as a kml file. 
%
%   This function receives a cell array of character arrays, each
%   containing a portion of the kml file to be written. The function
%   concatinates these arrays and then appends the array with a header and
%   footer to make the output a valid kml file. The returned char array can
%   then be written to a kml file. 

%INPUTS:
%flnm_list      nx1 cell array of char arrays containing portions of kml
%               data that will be concatinated
%docname        char array of the document name that will be written in the
%               kml file

%OUTPUTS:
%a char array that contains all the text of the combined kml files that can
%then be writted to a .kml file. 


header = sprintf(['<?xml version="1.0" encoding="utf-8"?> \r\n'...
                    '<kml xmlns="http://www.opengis.net/kml/2.2"> \r\n' ...
                    '\t<Document>\r\n'...
                    '\t<name>',docname,'</name>\r\n']);
      
tailer = sprintf(['\t</Document>\r\n'...
                  '</kml>']);
              
%The number of files to read
if iscell(charinputs)
    num_groups = length(charinputs);
else
    num_groups = 1;
    charinputs = {charinputs};% Change the form so we don't get a error with curly braces later
end
    

%Find the place in the first file where the document termination string
%occurs
text_hold = '';
for i = 1:num_groups
    text_i_select = charinputs{i};
    text_hold = [text_hold,text_i_select];%Add it to the compiled text string. 
end

%Merge the first file with all the others and stich on the document
%termination text to make it a valid kml. 
text_combined = [header,text_hold,tailer ];

end


      