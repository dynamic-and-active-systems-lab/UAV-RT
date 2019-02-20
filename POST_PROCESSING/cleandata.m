    function [data_out] = cleandata(Fs,pulse_dur,filt_freq,band,data_in)
%CLEANDATA bandpass filters input complex radio data and then finds the
%maximum of that filtered data at a sampling frequency equal to that of
% 1/pulse_duration

%   Radio data if fed to this function. The data is bandpassed filted with
%   a Chebyshev type II filter with a bandwidth and center frequency
%   defined in the input to the function. If the center frequency is less
%   than half the bandwidth of the filter, the low bound of the filter is
%   set to 0 Hz. The filtered data is then windowed at 2x size of the pulse
%   duration and the maximum is found in that window. A 50% overlap is used
%   when doing the windowed scan to be sure to capture an entire pulse in 
%   one of the windows. The output frequency is then 1/pulse_durs. 

%Conversions b/t Python and Matlab:     https://www.mathworks.com/help/compiler_sdk/python/pass-data-to-matlab-from-python.html

%INPUTS:
%TYPES        Python           Matlab          DESCRIPTION
%-------------------------------------------------------------------------
%Fs        -- int           --> int32/int64 -- sampling frequency of data_in(max is 2^32 = 4,294,967,295 Hz based on requirements for Python - Matlab type conversions. Int32 for windows. Int 64 for Linux and Mac)
%pulse_dur -- float         --> dbl         -- time of the pulse duration (s)
%filt_freq -- float         --> dbl         -- frequency center frequency of the bandpass filter(Hz)
%band      -- float         --> dbl         -- bandwidth of bandpass (Hz)
%filt_freq -- float         --> dbl         -- frequency targeted for pulse detection (Hz)
%data_in   -- matlab.double (nx1) --> numeric array  -- vector of complex radio data (in complex form with I being the real part and Q being the imaginary part in its complex representation)

%To generate the appropriate complex Python array is the is_complex option:
%a = matlab.double([[1+2j,2,3], [6,7,8]],is_complex=True)

%Python can only handle one output, so we pack the outputs into a cell
%array that Python reads in as a list. We can then unpack the list within
%Python.

%OUTPUT:
%TYPES        Matlab            Python                 DESCRIPTION
%-------------------------------------------------------------------------
%data_out   -- cell array (1,2) --> list   -- contains the time and data output vectors
%  |  
%  |
%  --> t_out    -- numeric array --> matlab.double (nx1) -- vector of times for data_out
%  |  
%  |
%  --> data_out -- numeric array --> matlab.double (nx1) -- vector of complex radio data


%% Format inputs
%This logic is only needed because sometimes we want to use this natively
%in Matlab and don't need the cell2mat command that is required when
%receiving a list from Python
Fs = double(Fs);
pulse_dur = double(pulse_dur);
filt_freq = double(filt_freq);
band = double(band);

if iscell(data_in)
    SDR_raw = cell2mat(data_in);
    CELL_FLG = 1;
else
    SDR_raw = data_in;
    CELL_FLG = 0;
end

%%

f_c = filt_freq;
half_band = band;  %we'll need this later when defining the filter


t_raw = 1/Fs*(0:1:length(SDR_raw)-1)+1/Fs; %create a time vector at raw sample rate. First sample is at 1/Fs
pulse_n = pulse_dur*Fs;             %the expected number of samples per pulse
T = length(SDR_raw)/Fs;             %the total time of the data

n_window = pulse_n*2;               %window of the sprectrogram is 2x the pulse duration
n_overlap = pulse_n;                %overlap of 50%


%% Filter the SDR data

low_freq = max([1,f_c-half_band]); %lower bound of bandpass. Max is used in case center frequency is less than 1/2 band away from 0 Hz. I use 1 Hz here not 0 Hz, because the filter can't accept 0 Hz
upp_freq = f_c+half_band;          %upper bound of bandpass

[b,a] = cheby2(3,20,[low_freq upp_freq]*(2/Fs),'bandpass'); %NOTE: if you want to do an 4th order or high filter the tf approach has numerical round off problems and you'll need to use the ss approach comment out below
%[b,a] = cheby2(3,20,[low_freq upp_freq]*2/Fs,'bandpass'); %NOTE: if you want to do an 4th order or high filter the tf approach has numerical round off problems and you'll need to use the ss approach comment out below
%[A,B,C,D] = cheby2(3,20,[f_c-half_band f_c+half_band]*2/Fs,'bandpass'); % I don't use the approach currently because 'lsim' is about 20 times slower than 'filter'
%SDR_filt_ss = lsim(ss(A,B,C,D),SDR_raw,t_raw);
SDR_filt = filter(b,a,SDR_raw);
    
%Sometimes we have datasets that we have stiched together with zeros. After
%filtering numerical errors create values that are non zero, but very low,
%like 1e-30 to 1e-300. The pulse finding algorithm can see variations in
%these errors as pulses, where there were none, so we just rezero the data
%here if it is really low. 
SDR_filt(abs(SDR_filt)<1e-6) = 0; 

%Process the filtered data to develop a metric of the pulse amplitude over
%each window. The window is as wide as the expected pulse duration.
    tick = 1;
    inds_for_windows = 1:n_overlap:length(SDR_raw)-n_window;
for i = inds_for_windows
    freq_indices = round(f_c/Fs*n_window) + 1;   
    data(tick) = max(SDR_filt(i:i+n_window-1));    %Using the max seems to work better
    %data(tick) = mean(SDR_filt(i:i+n_window-1));    %
    tick = tick+1;
end

%The g subscript is for the data afer resampling down to the width of the
%pulse repetition rate. 
number_of_windows = length(inds_for_windows);


Fs_g = 1/((n_window-n_overlap)/Fs);

t_g = 1/Fs_g*((0:1:length(data)-1))+(n_window/Fs*1/2);%shift by 1/2 the sampling period because the first measurement is at the center of the first window

%Helpful for debugging
    % figure
    % plot(t_raw,SDR_raw);hold on;
    % plot(t_raw,SDR_filt)
    % plot(t_g,data)


%% Format outputs

% if CELL_FLG == 1
%     data_out = mat2cell(double(data));
%     t_out    = mat2cell(double(t_g));
% else
%    data_out = data;
%    t_out = t_g;
%end

%data_out{1} = t_g;
%data_out{2} = data;

data_out{1,1} = double(t_g);
data_out{2,1} = double(data);



