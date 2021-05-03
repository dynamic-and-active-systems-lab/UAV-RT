function [OUTSTRUCT] = configread(filenamestr)
%CONFIGREAD reads and parses a UAV-RT flight configuration file and return
%a structure with the different fields specified in the config file. 
%
%
%%Inputs: 
%filenamestr    a string of the filename of the config file to parse
%
%Outputs:
%OUTSTRUC    a configuration structure including the following fields
%	-TS      Matlab dateime of the config file timestamp            
%	-TAG     1xn sturcture, where n is the number of tags specified in the
%               config file
%       -FT     scalar of the tag carrier frequency in MHz
%       -FRS    scalar sampling frequency of used by the radio (Hz)
%       -PUL_D  pulse duration of the tag in (s)
%       -PUL_R  pulse repetition period (s)
%   -SDRSYS  1x1 structure of SDR system configuration fields
%       -FVEH   scalar of the UAV telemetry sampling rate (Hz) by companion
%               computer
%       -RFG    scalar RF Gain in dB
%       -IFG    scalar IF Gain in dB
%       -BBG    scalar BB Gain in dB
%
%
%Author: Michael Shafer
%Date: 2019-05-30
%
%%
%Read the file
STRING_IN = fileread(filenamestr);
%Find all the new line character locations
NL_inds = regexp(STRING_IN,'\r\n');
%Find all the TAG character locations
TAG_inds = regexp(STRING_IN,'//TAG\r\n');
numoftags = length(TAG_inds);

%Here are all the fields we'll need to check. 
%Update this list later if you add mode fields to the config files 
SDRSYS_FIELDS = {'FVEH','RFG','IFG','BBG'};
TAG_FIELDS = {'FT','FRS','PUL_D','PUL_R'};

%Set up the time stamp
[~, e_ind] = regexp(STRING_IN,'TIME\t');
datestringin = STRING_IN(e_ind+1:-1+NL_inds(find(NL_inds>e_ind,1,'first')));
OUTSTRUCT.TS = datetime(datestringin,'InputFormat','dd-MMM-yyyy HH:mm:ss')';

%Set up the SDRSYS Fields
for i =1:length(SDRSYS_FIELDS) %Looking for all SDRSYS fields
    curr_field = SDRSYS_FIELDS{i};%The current field
    searchstr = [curr_field,'\t']; %The string of the field we are looking for
    [~, e_ind] = regexp(STRING_IN,searchstr); %Find all the places where this field is listed
    OUTSTRUCT.SDRSYS.(curr_field) = STRING_IN(e_ind+1:-1+NL_inds(find(NL_inds>e_ind,1,'first')));
end

%Set up the TAG Fields
for  i =1:length(TAG_FIELDS)  %Looking for all tag fields
    curr_field = TAG_FIELDS{i}; %The current field
    searchstr = [curr_field,'\t']; %The string of the field we are looking for
    [~, e_ind] = regexp(STRING_IN,searchstr); %Find all the places where this field is listed
    for j = 1:numoftags  %Now pull it out for each of the different tags
        OUTSTRUCT.TAG(j).(curr_field) = STRING_IN(e_ind(j)+1:-1+NL_inds(find(NL_inds>e_ind(j),1,'first')));
    end
end
%Set up the NOTES
[~,e_ind] = regexp(STRING_IN,'//NOTES\r\n');
OUTSTRUCT.NOTES = STRING_IN(e_ind+1:end); %Everything beyond //NOTES\n is in the notes

end