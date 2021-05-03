function [lossprcnt,cmdout] = isipup(IPaddress,n_pings,timeout,printechooption)
%isipup Pings the specified IP address to check if it is on the network
%   This program uses the system command to pin an IP address. n_pings
%   pings are made and the % of packet losses are reported as the quality
%   output. Timeout in seconds can be specified by the timeout setting. If
%   printechooption is set to true, the -echo command will be issued to
%   print the system terminal output the the Matlab command window. The
%   target IP address (IPaddress) variable shoudl be entered as a string
%
%Inputs:
%   IPaddress   The target IP address ex/ 10.0.1.5'
%   n_pings     Number of echo requests to send
%   timeout     Ping timeout of each echo request
%   printechooption  true or false (1,0) to indicate if cmd out should be
%                    printed to matlab terminal
%
%Outputs:
%   lossprcnt   The percentage of lost packets
%   cmdout      Raw output of the command window


%https://www.lifewire.com/ping-command-2618099
%http://www.tutorialspoint.com/unix_commands/ping.htm

if ismac
command = ['ping',...
            ' -c ', num2str(n_pings),...
            ' -t ', num2str(timeout*1000),...
            ' ',IPaddress];
elseif ispc
   command = ['ping',...
            ' -n ', num2str(n_pings),...
            ' -w ', num2str(timeout*1000),...
            ' ',IPaddress];
 
end 
if printechooption
    [~,cmdout] = system(command,'-echo');  
else
    [~,cmdout] = system(command);  
end

if ismac
    [startIndex,endIndex] = regexp(cmdout,',\s\w*.\w%'); %comma space number.number %
     lossprcnt = str2double(cmdout([startIndex+2:endIndex-1]));
else
    [startIndex,endIndex] = regexp(cmdout,'(\w*%\sloss');% left parentheses number % space loss
     lossprcnt = str2double(cmdout([startIndex+1:endIndex-6]));

end

