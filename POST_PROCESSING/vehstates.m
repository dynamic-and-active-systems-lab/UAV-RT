function [data_out] = vehstates(query_time_in,vehicle_states_in,waypt_time_in)
%VEHSTATES find the state of the vehilce at given instants in time
%   This function takes in a list of times and then finds the associated
%   vehicle states at those times. The first column of the state matrix
%   must be a list of times. It also determines a which waypoint the
%   vehicle was at the query times. If waypt_time is empty, the resulting
%   waypoint list is all NaN. 

%INPUTS:
%TYPES            Python           Matlab          DESCRIPTION
%-------------------------------------------------------------------------
%query_time_in -- matlab.double (nx1) --> numeric arry -- array of times where we want the vehicle state
%vehicle_states_in -- matlab.double (pxk) --> numeric arry -- array of vehicle states in time. Each k 
%                                                             column is a state. Each row is the total 
%                                                             state at a given instant in time. The first 
%                                                             column must be the time when the vehicle 
%                                                             state was assessed. 
%waypt_time_in -- matlab.double (mx2) --> numeric arry -- an array that has the times when the vehicle was
%                                                         at each waypoint. The row number corresponds to
%                                                         the waypoint number. THe first column is the
%                                                         time when the waypoint was arrived at and
%                                                         the second column is when that waypoint was
%                                                         departed. 
%                                                         
%                                                        
%                                                         

%ex/ %a = matlab.double([[1,2,3], [6,7,8]])

 

%Python can only handle one output, so we pack the outputs into a cell
%array that Python reads in as a list. We can then unpack the list within
%Python.

%OUTPUT:
%TYPES        Matlab            Python                 DESCRIPTION
%-------------------------------------------------------------------------
%data_out   -- cell array (1,2) --> list   -- contains the time and data output vectors
%  |  
%  |
%  --> query_states_out -- numeric array (nxk) --> matlab.double -- vehicle states at each n time to be assessed. This matrix includes the
%                                                                   query times in the first column so that number of columns of the vehicle
%                                                                   state matrix and the output are the same
%  |  
%  |
%  --> pulse_waypt_num_out  -- numeric array (nx1) --> matlab.double -- a list of which waypoint the pulse was received at. NaN values if not at a waypoint. 




%% Format inputs
%This logic is only needed because sometimes we want to use this natively
%in Matlab and don't need the cell2mat command that is required when
%receiving a list from Python

if iscell(vehicle_states_in)
    vehicle_states = cell2mat(vehicle_states_in);
else
    vehicle_states = vehicle_states_in;
end


if iscell(query_time_in)
    query_times = cell2mat(query_time_in);
else
    query_times = query_time_in;
end


if iscell(waypt_time_in)
    waypt_time = cell2mat(waypt_time_in);
else
    waypt_time = waypt_time_in;
end


%% Intperpolate vehilce states to pulse times and find average states

query_states = interp1(vehicle_states(:,1),vehicle_states(:,1:end),query_times); %INCLUDE THE TIME COLUMN
%avg_states = mean(vehicle_states(:,2:end));

%% FLIGHT DATA AND RADIO DATA FUSION

%This block of code determines at which waypoint each pulse was received.
%It looks at the time of each receive pulse and determines if it was in
%between time stamps when the vehicle was at a way point. If it was not at
%a waypoint, it puts NaN in the pulse_waypt_num vector. 

%pulse_waypt_num is a vector that is as long as the number of detected
%pulses and has a listing of the waypoint at which the pulse was detected

pulse_waypt_num = NaN*ones(size(query_times)); %Preallocate vector. 

if ~isempty(waypt_time)%If the request included times when the system was at a waypoint, then create a list of at which waypoint each pulse was received. 
    for i = 1:length(query_times)
        waypt_num_logic = (query_times(i)>waypt_time(:,1))&(query_times(i)<waypt_time(:,2));
        if max(waypt_num_logic)==1
            pulse_waypt_num(i) = find(waypt_num_logic==1,1,'first');
        else
            pulse_waypt_num(i) = NaN;
        end
    end
end


%% Format outputs to ensure compatibility with Python
query_states_out = double(query_states);
pulse_waypt_num_out =double(pulse_waypt_num);

data_out{1} = query_states_out;
data_out{2} = pulse_waypt_num_out;

end

