function [query_states_out,pulse_waypt_num] = vehstates(query_times,vehicle_states_in,waypt_time)
%function [data_out] = vehstates(query_time_in,vehicle_states_in,waypt_time_in)
%VEHSTATES find the state of the vehilce at given instants in time
%   This function takes in a list of times and then finds the associated
%   vehicle states at those times. The first column of the state matrix
%   must be a list of times. It also determines a which waypoint the
%   vehicle was at the query times. If waypt_time is empty, the resulting
%   waypoint list is all NaN. 

%INPUTS:
%TYPES               Matlab          DESCRIPTION
%-------------------------------------------------------------------------
%query_time     --> numeric arry -- array of times where we want the vehicle state
%vehicle_states_in --> numeric arry -- array of vehicle states in time. Each k 
%                          or          column is a state. Each row is the total 
%                       sturcture      state at a given instant in time. The first 
%                                      column must be the time when the vehicle 
%                                      state was assessed, or if passed as
%                                      a structure, the timestamps must be
%                                      labeled 'time' ex.
%                                      vehicle_states_in.time = [1 2..]
%
%waypt_time    --> numeric arry -- an array that has the times when the vehicle was
%                                     at each waypoint. The row number corresponds to
%                                     the waypoint number. The first column is the
%                                     time when the waypoint was arrived at and
%                                     the second column is when that waypoint was
%                                     departed. 

%OUTPUT:
%TYPES        Matlab                       DESCRIPTION
%-------------------------------------------------------------------------
%query_states_out -- numeric array (nxk) -- vehicle states at each n time to be assessed. This matrix includes the
%                          or               query times in the first column so that number of columns of the vehicle
%                     structure             state matrix and the output are
%                                           the same. If a structure was
%                                           provided, the output structure
%                                           names are retained.
%pulse_waypt_num_out  -- numeric array (nx1) -- a list of which waypoint the pulse was received at. NaN values if not at a waypoint. 




%% Format inputs
%This logic is only needed because sometimes we want to use this natively
%in Matlab and don't need the cell2mat command that is required when
%receiving a list from Python

%If input was a structure, make the vehicle states a matrix to make the
%interpolation code cleaner. We'll repack into an output structure after
%interpolating.

if isstruct(vehicle_states_in)
    vehicle_states(:,1) = vehicle_states_in.time;
    fn      = fieldnames(vehicle_states_in);
    num_fn  = numel(fn);
    tick = 2; %Tick specifies the column placement of the state vector in
    %the state matrix. We start at two because the first column
    %is always used for time.
    column_placement = zeros(1,num_fn);%Will be as long as the number of fieldnames
    for i = 1:num_fn
        if ~(strcmp(fn{i},'time')) %We'll skip over time becasuse it was already placed
            vehicle_states(:,tick) = vehicle_states_in.(fn{i});
            column_placement(i) = tick; %Keep a record of in which
            %column the ith field was
            %placed so that we can repack a
            %new output structure later.
            tick = tick+1; %Only update the column placement if we
            %used the last field. We we encounter the
            %time field, we update i, but not tick
        else %If we encounter the time field, we record that it's
            %placement was in the first column
            column_placement(i) = 1;
        end
    end
else
    vehicle_states = vehicle_states_in;
end

%Vestigial code
% if iscell(vehicle_states_in)
%     vehicle_states = cell2mat(vehicle_states_in);
% else
%     vehicle_states = vehicle_states_in;
% end
% 
% 
% if iscell(query_time_in)
%     query_times = cell2mat(query_time_in);
% else
%     query_times = query_time_in;
% end
% 
% 
% if iscell(waypt_time_in)
%     waypt_time = cell2mat(waypt_time_in);
% else
%     waypt_time = waypt_time_in;
% end


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


%% Format outputs to be similar to inputs. 
%If a sturcture was provided as an input, then provide a similar
%strucuture, but at the queried times. If a matrix was provided, provide
%the matrix output.
%
if isstruct(vehicle_states_in)
    for i = 1:num_fn
        query_states_out.(fn{i}) = query_states(:,column_placement(i));
    end
else
        query_states_out = query_states;
end


end

