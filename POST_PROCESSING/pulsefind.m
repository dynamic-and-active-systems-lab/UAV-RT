function [data_out] = pulsefind(Fs,pulse_dur,pulse_rep,data_in)
%PULSEFIND generates a list of pulse amplitude and times based on inputed
%radio data.
%   This function takes in radio data at a given sampling frequency with an
%   expected pulse repetition rate. The algorithm first takes the absolute
%   value of the data in case complex data was passed to it. A moving
%   average of the signal is then taken so that the peaks can be compared
%   to moving average. A threshold is then used to determine a set of
%   candidate peaks for consideration as pulses. 
%
%   A square wave with a frequency of 1/pulse_rep is then cross correlated 
%   with data to eliminated those peaks between pulses. The square wave 
%   cycles between 0 and 1 with a 10% duty cycle. In theory, the duty cycle
%   should only be pulse_dur/pulse_rep*100, but we don't know the pulse_rep
%   exactly so we use 10%. After the appropraite shift on the square wave 
%   is determined using the cross correlation, the square wave and the data
%   is then multiplied to zero out false peaks between currently believed 
%   true peaks. The cross correlation is done using three cycles of the 
%   square wave and the equivalen chunk of the data. Because we don't know 
%   the pulse_rep exactly, if you attempt to cross correlate the entire
%   dataset, the slight difference in frequency will result in periodic
%   fading of detected pulses. 

%   After the cross correlation eliminated obvious false positives between
%   peaks, we are left with a dataset of peaks that are spaced by
%   approximately pulse_rep aparte in time. There are however adjacent
%   peaks that were allowed through because of the 10% duty cycle of the
%   square wave. They my be spaced out by as little as pulse_dur. This
%   likely results because of the windowing used when cleanind the original
%   data. That windowing scheme uses a 2xpulse_dur window with 50% overlap.
%   The result is that a single pulse could show up 3 times in the log, but
%   only one of those detection would contain the entire pulse. To deal
%   with this, this algorithm looks through the list of potential pulses
%   and if there are some that are spaced in time less tha 20% of the pulse
%   repetition rate, they are flagged and only that pulse within that group
%   with the maximum amplitude is logged as the real pulse. After these
%   adjacent pulses are consolidated/eliminated the pulse list and times
%   when those pulses occured are output.


%INPUTS:
%TYPES        Python           Matlab          DESCRIPTION
%-------------------------------------------------------------------------
%Fs        -- int           --> int32/int64 -- sampling frequency of data_in(max is 2^32 = 4,294,967,295 Hz based on requirements for Python - Matlab type conversions. Int32 for windows. Int 64 for Linux and Mac)
%pulse_dur -- float         --> dbl         -- time of the pulse duration (s)
%pulse_rep -- float        	--> dbl         -- repetition rate of pulses (s)
%data_in   -- matlab.double (nx1) --> numeric array  -- vector of complex radio data

%To generate the appropriate complex Python array is the is_complex option:
%a = matlab.double([[1+2j,2,3], [6,7,8]],is_complex=True)



%OUTPUT:
%TYPES        Matlab            Python                 DESCRIPTION
%-------------------------------------------------------------------------
%data_out   -- cell array (1,2) --> list   -- contains the time and data output vectors
%  |  
%  |
%  --> pulse_time_out -- numeric array --> matlab.double (nx1) -- vector of times when pulses were detected
%  |  
%  |
%  --> pulse_list_out -- numeric array --> matlab.double (nx1) -- vector of pulse amplitudes



%OUTPUT:
%TYPES        Matlab            Python          DESCRIPTION
%-------------------------------------------------------------------------
%
%


Fs = double(Fs);
pulse_dur = double(pulse_dur);
pulse_rep = double(pulse_rep);


if iscell(data_in)
    data_in_2 = cell2mat(data_in);
else
    data_in_2 = data_in;
end


t = 1/Fs * (0:1:length(data_in_2)-1)+1/Fs; %First point is at 1/Fs


data_abs = abs(data_in_2);


%% CREATE A MOVING MEAN OF THE DATA TO USE AS A BASELINE FOR THE THRESHOLD DETECTOR
%use moving mean over 1 pulse rep rate
n_pulse_rep = round(pulse_rep*Fs);%this is the number of samples at the windowed rate that should containt a single pulse based on the rep rate of the pulses
n_pulse_dur = round(pulse_dur*Fs);
movemean_data = movmean(data_abs,n_pulse_rep-5*n_pulse_dur);%moving mean of the data
%I want the moving mean to only ever include 1 pulse, so I use the pulse
%repetition rate as the window size on the moving mean, but shorten it by
%5x the pulse duration so that when the second pulse hits, the first is no
%longer included in the current moving mean value. 

% figure;
% plot(t,data_abs); hold on;
% plot(t,movemean_data); hold on;


%% THRESHOLD TO DETECT PULSES OVER 1.5X THE MOVING MEAN
data_abs_thresh = data_abs;
data_abs_thresh(data_abs<1.5*movemean_data)=0;
% plot(t,data_abs_thresh,'--'); hold on;

%% SLIDING WINDOW CORRELATOR TO REJECT FALSE POSITIVES IN THE TIME BETWEEN PULSES
tick = 1;
pulse_sunset = 3;   %How many pulses in the past do you want to compare against in the sliding correlator? 
n_window = round(Fs*pulse_rep*pulse_sunset);    %look back for pulse_sunset pulse periods
n_pulse_rep = round(Fs*pulse_rep);              %width of a single pulse period
t_window = pulse_rep*pulse_sunset;                %time of the window we are considering (time of the pulse_sunset window)

for i = n_window+1:n_pulse_rep:length(data_abs_thresh)  %we jump forward by n_pulse_rep each time and then look back n_window
    data_abs_thresh_wind = data_abs_thresh(i-n_window:i); %this is the block of data we want to consider
    %This next line creates a square wave with 0-1 amplitude that has a
    %period of the pulse_rep and a duty cycle of 10% so that we can be off
    %by +/-5% of the actual pulse rep rate and still catch the pulses in
    %the sliding correlator. 
    y = 1/2*(1+square(2*pi*(1/pulse_rep)*t(i-n_window:i),10));%+/- 5% on the pulse
    [xcorvals, thelags]=xcorr(data_abs_thresh_wind,y);  %Do the sliding window correlation
    shift_inds = thelags(find(xcorvals==max(xcorvals),1,'first'));  %Find the index shifting that has the highest correlation
    y_shift = circshift(y,shift_inds);  %Now shift the square wave by the shift index that maximized correlation
    data_abs_thresh_xcor(i-n_window:i) = data_abs_thresh_wind.*y_shift;   %Multiply the data and the shifted square wave to silence spurrious signals between actual pulses
    y_shift_log(i-n_window:i) = y_shift; %Take the current processed block and add it to the record
    tick = tick+1;
end

if length(data_abs_thresh_xcor)<length(t) %in most cases the previous block will have a few less points than the original vector, so we just set the last few of the xcor vector equal to the original
    data_abs_thresh_xcor(end+1:length(t)) =  data_abs_thresh(length(data_abs_thresh_xcor)+1:length(t));
    y_shift_log(end+1:length(t)) =  zeros(1,length(t)-length(y_shift_log));
end

%plot(t,data_abs_thresh_xcor,':c'); hold on;
%plot(t,y_shift_log*0.25,':k'); hold on;

%Now create a list of pulses and their times. 
pulse_list_1 = data_abs_thresh_xcor(data_abs_thresh_xcor~=0);
pulse_time_1 = t(data_abs_thresh_xcor~=0);

%% ADJACENT PUSLE DETECTION
%THIS BLOCK OF CODE IS USED TO ELIMINATE PULSES THAT WERE IMMEDIATELY
%ADJACENT TO ONE ANOTHER. IT BASICALLY COMPARES THE TIME STAMPS OF ADJACENT
%PUSLES AND IF THEY ARE TOO CLOSE TOGETHER IN TIME, IT TAKES THE AVERGAGE OF
%THE TWO (or more). THE RESULTING PULSE LIST AND TIME FOR THOSE PUSLES ARE SAVES AS A
%NEW LIST OF PULSES.


%I uses this matrix to think throught this algorithm
%themat = [pulse_list_1',pulse_time_1',[(diff(pulse_time_1)<0.2*pulse_rep)';0]];
dbl_cnt_logic = [(diff(pulse_time_1)<0.2*pulse_rep)';0];
tick = 1;
i = 1;
while i < length(pulse_list_1)   %We are running through the entire pulse detection list
    %disp(['i = ',num2str(i)])
    if dbl_cnt_logic(i) == 1    %First check to see if the dt between pulses was less than the threshold
        click = 1;              %Click keeps track of advancing indecies beyond i for looking forward to the pulses that had a dt between the current pulse and their time less than the threshold
        avg_list = [pulse_list_1(i),pulse_list_1(i+1)]; %Create a list of the ith and ith+1 pulse
        avg_time = [pulse_time_1(i),pulse_time_1(i+1)]; %and their time
        while dbl_cnt_logic(i+click)==1         %Now check if the i+1 pulse was had a dt to the ith+2 pulse less than the threshold. If so, the ith ith+1 and ith+2 should all be averaged together then.
     %       disp(['click = ',num2str(click)])
            avg_list = [avg_list, pulse_list_1(i+click+1)]; %Add that next pulse to the list to be averaged
            avg_time = [avg_time, pulse_time_1(i+click+1)];
            click = click + 1;                  %Advanced the click to check the ith+3 (and so on) if necessary
        end
        %pulse_list_2(tick) = mean(avg_list);    %Take the average of the list as the real value
        %pulse_time_2(tick) = mean(avg_time);
        [pulse_list_2(tick), max_ind]  = max(avg_list);    %Take the max of the list as the real value - We used a 2x pulse length window with 50% overlap, so one of them will have the entire pulse
        pulse_time_2(tick) = avg_time(max_ind);             %Use the associated time of the max pulse
        
        i = i+2+click-1;                        %Advance the i counter to skip over all the pulses that were just averaged
    else
        pulse_list_2(tick) = pulse_list_1(i);   %If the dt was less than the threshold, then just take the current pulse as the pulse
        pulse_time_2(tick) = pulse_time_1(i);
        i = i+1;
    end
    tick = tick+1;                              %Advance the counter of the cleaned up pulse list
end


pulse_time_out = double(pulse_time_2);
pulse_list_out = double(pulse_list_2);

%data_out{1} = pulse_time_out;
%data_out{2} = pulse_list_out;

data_out{1} = pulse_time_out;
data_out{2} = pulse_list_out;

end

