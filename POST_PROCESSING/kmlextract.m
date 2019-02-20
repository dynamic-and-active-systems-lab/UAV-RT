function [extract] = kmlextract(flnm)
%KMLEXTRACT pulls and returns just the contents of a kml file, not include
%text before the first </name> string, nor after the </Document> string. 
%   This function receives an single kml filename extension. From that
%   files it pulls the contents after the first </name> and before the last
%   <\Document> as a string

%INPUTS:
%flnm           Char array of the .kml filename. Ex: 'data.kml'

%OUTPUTS:
%extract        char array extracted from original kml


%Find the place in the first file where the document termination string
%occurs
    texti = fileread(flnm);
    head_1_ind = min(strfind(texti,'</name>')+7);
    tail_1_ind = strfind(texti,'</Document>')-1;
    extract = texti(head_1_ind:tail_1_ind); %Pull out the content
end

