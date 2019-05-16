
function [out] = radioprep(data_filename,varargin)
%RADIO_PREP reads in a IQ data .dat file produced by GNU radio and finds
%the frequency of pulses within that radio data if they are below 1kHz
%   This function  first reads in the IQ data, then converts it to a
%   complex representation. After the conversion, a subset of the data (if
%   timeselect is used) is Fourier transformed to find the double-sided
%   amplitude. This is then smoothed with a moving average over frequency
%   with 1Hz windows. The peak in the result is then found between -f_max
%   and f_max and reported as the pulse frequency. The spectrum can be 
%   plotted with options, and a spectrogram can be plotted with options. 
%
%Required Inputs:
%   data_filename       char array filename or path to .dat radio data log of
%                       IQ radio data saved by GNU radio
%
%Options (don't need to spell out the name if value comes right after filename)
%   'fs'                sampling rate of radio data in Hz. Default is 48000
%                       if not entered
%   'pulse_dur'         velocity (m/s) threshold below which you can be 
%                       considered at a waypoint. Default is 0.5 m/s if
%                       not entered.
%   'timeselect'        [t_start t_end] start time and end time that you
%                       want to consider, measured from the start of the 
%                       record. This doesn't change the length of the
%                       output, only the region overwhich you want to use
%                       to conduct the FFT looking for peaks
%Parameters (need to spell out the name, then enter the value)
%   'plot'              Control of plotting of results. Enter 'fft' to plot
%                       FFT of entire dataset. Enter 'spectro' to see the
%                       spectrogram. Enter 'both' for both plots. Default
%                       is no plotting. 
%   'filterband'        band of planned filter. Default is 160 Hz. Used
%                       only to plot the band in the results, if plotting 
%                       is turned on. 
%   'fmax'              The maximum frequency from base band where tag
%                       frerquency may be expected. Default is 1 kHz
%   'axis_handles'      Handle of axis where data should be plotted. If
%                       empty, function will create a new figure. If
%                       plotting both, the FFT should be the first handle
%                       and spectrogram axes should be the second. ie. 
%                        axis_handles = [ax_fft, ax_spectro];


%                         
%Outputs:
%   'out'                 a cell array with the following elements
%                       out{1}  f   The found pulse frequency
%                       out{2}  SDR_raw raw IQ data in complex form, same
%                               length as the input data
%                       out{3}  Frequency vector from the FFT
%                       result
%                       out{4}  Singled-sided amplitude results from FFT
%
%


%% Parse the input variables
%Modeled from 
% https://www.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
p = inputParser;

default_fs = 48000;
default_pulse_dur = 0.02; %20 ms default pulse duration
default_filt_band = 160; %160 Hz filter band
default_plot = 'none';%Don't plot by default
default_t_bounds = 0; %Set to zero so we can use the bounds of the time vector later if nothing is entered. 
default_f_max = 1000; %Set to zero so we can use the bounds of the time vector later if nothing is entered. 
default_axis = [];   %Create empty variable for figure handle. 

checkplot = @(x) any(validatestring(x,{'fft','spectro','both','none'}));
allhandle = @(x) all(ishandle(x))||isempty(x); %Check all array elements are handles or is empty for default case

%Required inputs
addRequired(p,'file_location',@ischar);
%Not required and doesn't need to spell out .....,[10 30]) if used after
%the filename, fs, pulse_dur
addOptional(p,'fs',default_fs,@isnumeric);
addOptional(p,'pulse_dur',default_pulse_dur,@isnumeric);
addOptional(p,'timeselect',default_t_bounds,@isnumeric);
addOptional(p,'fmax',default_f_max,@isnumeric);
%Not required but must be spelled out like .....'timeselect',[10 30]) if used
addParameter(p,'plot',default_plot,checkplot);
addParameter(p,'filterband',default_filt_band,@isnumeric);
addParameter(p,'axis_handles',default_axis,allhandle);%If array, all must be a handle


parse(p,data_filename,varargin{:})

%Useful for parser debugging
% if ~isempty(fieldnames(p.Unmatched))
%    disp('Extra inputs:')
%    disp(p.Unmatched)
% end
% if ~isempty(p.UsingDefaults)
%    disp('Using defaults: ')
%    disp(p.UsingDefaults)
% end

telemfile = p.Results.file_location;
Fs = p.Results.fs;
pulse_dur = p.Results.pulse_dur;
t_bounds = p.Results.timeselect;
filt_band = p.Results.filterband;
plot_type =  p.Results.plot;
f_max =  p.Results.fmax;
ax2plot =  p.Results.axis_handles;

%% Get plot axes ready
%Don't need to do anything if 
if strcmp(plot_type,'both')
    if isempty(ax2plot)%Create the 2 figures and get the axes handles for later
        ax2plot = gobjects(1,2);%Create the axes array for population later
        figure;
        ax2plot(1) = gca;
        figure;
        ax2plot(2) = gca;
    else
        if length(ax2plot)~=2
            error('Two axis handles needed if plotting both.')
        end
    end
elseif strcmp(plot_type,'fft')
    if isempty(ax2plot)%Create the 1 figures and get the axes handles for later
        ax2plot = gobjects(1,2);%Create the axes array for population later
        figure;
        ax2plot(1) = gca;
    else
        if length(ax2plot)~=1
            error('Only one axis handle needed.')
        else %If they passed only one handle, but it in the right place in the ax2plot array
            ax2plot_hold = gobjects(1,2);%Create the axes array for population later
            ax2plot_hold(1) = ax2plot;
            ax2plot = ax2plot_hold;
        end
    end
elseif strcmp(plot_type,'spectro')
    if isempty(ax2plot)%Create the 1 figures and get the axes handles for later
        ax2plot = gobjects(1,2);%Create the axes array for population later
        figure;
        ax2plot(1) = [];
        ax2plot(2) = gca;
    else
        if length(ax2plot)~=1
            error('Only one axis handle needed.')
        else
            ax2plot_hold = gobjects(1,2);%Create the axes array for population later
            ax2plot_hold(2) = ax2plot;    
            ax2plot = ax2plot_hold;
        end
    end
    
end    

%% Begin computation

count=inf;
f = fopen (telemfile, 'rb'); 
d = fread (f, [2, count], 'float'); %read in img and real parts of file
fclose (f);
SDR_raw = (d(1,:) + d(2,:)*1i)'; %combine real and imag into one vector and transpose to make a column
clear d f;


%SDR_raw 
%This is the raw SDR data with complex form with I being the real part and Q being the imaginary part in its complex representation
%https://www2.nau.edu/uavrt-p/wp-content/uploads/2017/06/I_Q_Data_Documentation_v1.pdf

%Fs = 48000;         %raw sampling rate (Hz)
%pulse_dur = 0.02;   %expected pulse duration (s)
%pulse_rep = 1.3;    %expected pulse repetition rate (s)

%t_raw = 1/Fs*(0:1:length(SDR_raw)-1); %create a time vector at raw sample rate
%T = length(SDR_raw)/Fs;             %the total time of the data

if t_bounds==0 %If user didn't specify a time bounds, use the entire dataset. 
    t_bounds = [0 length(SDR_raw)/Fs];
end

%Use Fs and start and end times to down select data to FFT. 
ind_start = floor(min(t_bounds)*Fs+1);%+1 because first element is at t = 0;
ind_end = floor(max(t_bounds)*Fs);
inds_select =  ind_start:1:ind_end;
if mod(length(inds_select),2) %Then it is odd
    inds_select = inds_select(1:end-1); %we want an even number of elements for future operations.
end

the_data = SDR_raw(inds_select);
%the_data = real(the_data)-1i*imag(the_data);

if strcmp(plot_type,'spectro')||strcmp(plot_type,'both') == 1
  pulse_n = pulse_dur*Fs;             %the expected number of samples per pulse
  %pulse_rep_n = pulse_rep*Fs;             %the expected number of samples per pulse
  %[s_spec,f_spec,t_spec] = spectrogram(the_data,pulse_n*2,pulse_n,1:10:10000,Fs);
  [s_spec,f_spec,t_spec] = spectrogram(the_data,pulse_n*2,pulse_n,-f_max:10:f_max,Fs);
  %[s_spec_2,f_spec_2,t_spec_2] = spectrogram(SDR_raw,pulse_n,0,1:10:10000,Fs);
  %[s_spec_3,f_spec_3,t_spec_3] = spectrogram(SDR_raw,pulse_rep_n*2,0,1:10:10000,Fs);
  % figure; surf(t_spec,f_spec,abs(s_spec),'Edgecolor','none');xlabel('Time(s)'); ylabel('Frequency (Hz)')
  surf(ax2plot(2),t_spec,f_spec,10*log10(abs(s_spec)),'Edgecolor','none');
  xlabel(ax2plot(2),'Time(s)'); 
  ylabel(ax2plot(2),'Frequency (Hz)')
  view(ax2plot(2),0,90)
end

%% SINGLE SIDED SPECTRUM
%  Y = fft(the_data);
%  L = length(the_data);
%  P2 = abs(Y/L);
%  P1 = P2(1:L/2+1);
%  P1(2:end-1) = 2*P1(2:end-1);
%  freqs = Fs*(0:(L/2))/L;
% figure; hold on
%  plot(freqs,P1,'b');
%  
%% DOUBLE SIDED SPECTRUM
 n = 2^nextpow2(size(the_data,1));
 L = length(the_data);
 Y = fft(the_data,n);
 P2 = abs(Y/L);
 P1 = [P2(n/2+2:end); P2(1:n/2+1)]; %This is just the two sided oriented. I'm using P1 because it is was what in the code later on from when I was only doign the signle sided spectrum
 freqs = Fs/2*linspace(-1,1,n);
 
 %figure;
 %plot(freqs_d,[P2_d(n/2+2:end); P2_d(1:n/2+1)],'--r')
  
 P1_smooth = movmean(P1,floor(length(freqs)/(freqs(end)-freqs(1))));
 freqs_1Hz = freqs(1):1:freqs(end);
 P1_1Hz = interp1(freqs,P1_smooth,freqs_1Hz);
 
 %f_max = 1000;%Maximum distance from center frequency to expected pulse frequency
 sub_range_select = (freqs_1Hz<f_max & freqs_1Hz>-f_max);
 freqs_1Hz_sub_range = freqs_1Hz(sub_range_select);
 [max_pow_val,ind_max_pow] = max(P1_1Hz(sub_range_select));%Find the maximum in the frequency range less than 1kHz- the tag should be within 1 kHz if the radio was tuned to its center frequency since with have 3 digits past the MHz decimal place - ie. 150.256 MHz - maxx it could be off is +/- 500 Hz 
 
 f = freqs_1Hz_sub_range(ind_max_pow); %TARGET FREQUENCY
 

 if strcmp(plot_type,'fft')||strcmp(plot_type,'both') == 1
     %figure
%      subplot(3,1,1)
%      t = ((0:length(SDR_raw)-1)/Fs)';
%      plot(t(inds_select),real(the_data));hold on
%      plot(t(inds_select),imag(the_data));     
%      title('IQ data in time domain')
%      xlabel('time (s)')
%      ylabel('IQ amplitude')
%      
%      subplot(3,1,2)
%      plot(freqs_1Hz,P1_1Hz)
%      title('Double-Sided Amplitude Spectrum of X(t)')
%      xlabel('f (Hz)')
%      ylabel('|P1(f)|')
% 	 set(gca,'ylim',[0 1.2*max_pow_val])
%      hold on
%      amp_range = get(gca,'ylim');hold on;
%      fill([f-filt_band/2,f-filt_band/2,f+filt_band/2,f+filt_band/2],[amp_range flip(amp_range)],[1 0.3 0.3],'Facealpha',0.1,'Edgecolor','none')
%      plot(freqs_1Hz_sub_range(ind_max_pow),max_pow_val,'r.','Markersize',10)
% 	 legend('Amplitude spectrum','Selected frequency')
% 
%      
%      subplot(3,1,3)
     plot(ax2plot(1),freqs_1Hz,P1_1Hz)
     title(ax2plot(1),['Double-Sided Amplitude Spectrum of X(t) Zoomed to \pm',num2str(-f_max),'Hz'])
     xlabel(ax2plot(1),'f (Hz)')
     ylabel(ax2plot(1),'|P1(f)|')
     set(ax2plot(1),'xlim',[-f_max f_max])
     set(ax2plot(1),'ylim',[0 1.2*max_pow_val])
     hold(ax2plot(1),'on')
     amp_range = get(ax2plot(1),'ylim');
     fill(ax2plot(1),[f-filt_band/2,f-filt_band/2,f+filt_band/2,f+filt_band/2],[amp_range flip(amp_range)],[1 0.3 0.3],'Facealpha',0.1,'Edgecolor','none')
     plot(ax2plot(1),freqs_1Hz_sub_range(ind_max_pow),max_pow_val,'r.','Markersize',10)
     legend(ax2plot(1),'Amplitude spectrum','Selected frequency')
 end
 
 out{1} = f;
 out{2} = SDR_raw;
 out{3} = freqs_1Hz;
 out{4} = P1_1Hz;
 
end

