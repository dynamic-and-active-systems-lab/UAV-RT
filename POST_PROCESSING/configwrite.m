function [] = configwrite(filenamestr,INSTRUCT,overwritecheck)
%CONFIGWRITE write UAV-RT flight configuration file.
%
%Inputs:
%filenamestr    a string of the filename of the config file to parse
%
%INTSTRUC    a configuration structure including the following fields
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
%overwritecheck     A string of either 'check' or 'overwrite'. This is used
%                   to automatically overwrite the existing file, if it
%                   already exists, or to check with the user for
%                   confirmation/new filename
%
%Note that the \r\n characters are used to specify new lines so that
%Windows systems can display the config files properly in Notepad.
%
%Author: Michael Shafer
%Date: 2019-05-30
%
%%
% filenamestr = 'config.uavrt';
% INSTRUCT.TS = datestr(now);
%
% INSTRUCT.SDRSYS.FVEH = 10;
% INSTRUCT.SDRSYS.RFG = 7;
% INSTRUCT.SDRSYS.IFG = 13;
% INSTRUCT.SDRSYS.BBG = 10;
%
% INSTRUCT.TAG(1).FT = 148.356;
% INSTRUCT.TAG(1).FRS = 48000;
% INSTRUCT.TAG(1).PUL_D = 0.02;
% INSTRUCT.TAG(1).PUL_R = 1.3;
%
% INSTRUCT.TAG(2).FT = 149.356;
% INSTRUCT.TAG(2).FRS = 52000;
% INSTRUCT.TAG(2).PUL_D = 0.05;
% INSTRUCT.TAG(2).PUL_R = 1.0;
%
% INSTRUCT.NOTES = 'These are notes provided about the flight.';

if ismac
    slash_char = '/';
elseif ispc
    slash_char = '\';
end

%First determine if a path is
slash_all_inds = strfind(filenamestr,slash_char);

if isempty(slash_all_inds)
    %The file is specified in the current dir. Check that dir for exisiting
    %filenames
    dest_dir = dir;%Current directory
else
    slash_ind = slash_all_inds(end);
    %The destination is in a different place
    dest_path = filenamestr(1:slash_ind-1);   %The path to the destination
    dest_name = filenamestr(slash_ind+1:end); %The file name
    dest_dir = dir(dest_path);            %Destination directory
end

%Check to see if any files exist in this destination with the same 
%filename. If so, check with the user to see if overwriting is okay.
if strcmp(overwritecheck,'check')%Only check with user is asked to
    for i = 1:length(dest_dir)
        if strcmp(dest_dir(i).name,dest_name)%Check to see if any
            msg = {['There already exists a configuration file in this',...
                ' folder with the name ',dest_name,'.'],...
                '',...
                'Do you want to overwrite or specify a new filename?'};
            answer = questdlg(msg,'Overwrite file?','Cancel','Overwrite','New Filename','New Filename');
            if strcmp(answer,'New Filename')
                default_new_name = [dest_name(1:end-6),'(2)',dest_name(end-5:end)];
                fileanswer = inputdlg('Specify filename','Filename',[1 40],{default_new_name});
                if ~isempty(fileanswer)
                dest_name = fileanswer{1};
                %Create a new name to open. 
                filenamestr = [dest_path,slash_char,dest_name];
                else%The user canceled
                    answer = 'Cancel';%Will stop the remainder of function
                end
            end
        end
    end
else
    answer = 'Overwrite';
end

if ~strcmp(answer,'Cancel')%Only do this if the user didn't cancel.
    
    fileID = fopen(filenamestr,'w');
%Write the creation timestamp field    
    fprintf(fileID,...
        ['//TS\r\n',...
        'TIME\t%s\r\n\r\n'],...
        INSTRUCT.TS);
    
%Write the SDRSYS fields    
    fprintf(fileID,...
        ['//SDRSYS\r\n',...
        'FVEH\t%2.1f\r\n',...
        'RFG\t%2.1f\r\n',...
        'IFG\t%2.1f\r\n',...
        'BBG\t%2.1f\r\n\r\n'],...
        INSTRUCT.SDRSYS.FVEH,...
        INSTRUCT.SDRSYS.RFG,...
        INSTRUCT.SDRSYS.IFG,...
        INSTRUCT.SDRSYS.BBG);

%Write the TAG fields
    for i = 1:length(INSTRUCT.TAG)
        fprintf(fileID,...
            ['//TAG\r\n',...
            'NUM\t%02u\r\n',...
            'FT\t%4.6f\r\n',...
            'FRS\t%-6u\r\n',...
            'PUL_D\t%2.3f\r\n',...
            'PUL_R\t%2.3f\r\n\r\n'],...
            i,...
            INSTRUCT.TAG(i).FT,...
            INSTRUCT.TAG(i).FRS,...
            INSTRUCT.TAG(i).PUL_D,...
            INSTRUCT.TAG(i).PUL_R);
    end
%Write the notes    
    fprintf(fileID,...
        ['//NOTES\r\n',...
        '%s'],...
        INSTRUCT.NOTES{1});
    
    fclose(fileID);
end
end
