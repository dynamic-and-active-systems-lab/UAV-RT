function [f_report,PRI,PD,PRI_std,PD_std] = tagcheck(radioin,Fs)
%TAGCHECK extracts the critical pulse characteristics from a tag.
%   This program is used to find the characteristics of a VHF tag.
%   Specifically, the carrier frequency difference from baseband, the pulse
%   repetition interval, and the pusle duration, all of which are necessary
%   for detection. The input radio data should be collected in a controlled
%   environment where there is a very high signal to noise ratio. The
%   algoithm used here assume very little noise. 
%
%   This function receives complex IQ radio data (radioin) at a given
%   sample rate (Fs). It then finds the frequency with maximum amplitude as
%   provided by an FFT. The program selects this frequency as the carrier
%   frequency (delta from baseband) as that of the tag. A third order
%   bandpass butterworth filter is then built around that frequency (low
%   pass if too close to base band). The data is filtered with a zero-phase
%   shift filter. The output is then shifted to baseband and lowpass
%   filtered again to eliminate higher fequency component introduced in the
%   frequency shifting. The amplitude of this signal is then thresholded
%   above 10% of the maximum in the entire waveform to indicate when pulses
%   are present. Any transitions to true (pulse active) that occur within 5
%   ms of a previousl transition to false (pusle inactive) are removed, as
%   this is an artifact of sinusoidal nature of the trailing edge of the
%   pulse logic wayforms. From this cleaned logic, transitions that took
%   more than 0.5 s are considered dwell times (inactive pulse).
%   Transitions that took less than 0.5 s, but more than 10ms, are
%   considered active pulse times. From this record, average and standard
%   deviations of pulse repetition interval and pulse durations are
%   calculated.
%
%INPUTS:
%   radioin     nx1 complex         Complex IQ input
%   Fs          1x1 scalar          Samplerate of radio in (Hz)
%
%OUTPUTS:
%   f           scalar              Delta frequency from baseband of
%                                   carrier frequency of pulse
%   PRI         scalar              Mean pulse repetition interval (s)
%   PRI_std     scalar              Standard deviation of pulse repetition 
%                                   interval (s)
%   PD          scalar              Mean pulse duration interval (s)
%   PD_std      scalar              Standard deviation of pulse duration 
%                                   (s)
%
%Author: Michael W. Shafer
%Date: 2019-08-12
%

%% FIND CARRIER

%FFT of input
n = 2^nextpow2(size(radioin,1));
L = length(radioin);
Y = fft(radioin,n);
P2 = abs(Y/L);
P1 = [P2(n/2+2:end); P2(1:n/2+1)]; %This is just the two sided oriented. I'm using P1 because it is was what in the code later on from when I was only doign the signle sided spectrum
freqs = Fs/2*linspace(-1,1,n);
%figure;plot(freqs,P1)

f_nyq = Fs/2; %Nyquist frequency

[max_pow_val,~] = max(P1(abs(freqs)<0.8*f_nyq));%Find the maximum in the frequency range below 90% of the nyquist
ind_max_pow = find(P1==max_pow_val,1,'first'); %Can't use the ind of the output of max from last line, because indicies are for the vector that is shortened to only consider freqs below 80% of nyquist
f_report = freqs(ind_max_pow); %%TARGET FREQUENCY
f = abs(f_report); %ABS TARGET FREQUENCY - needed for filtering


%% FILTER AND FEWUENCY SHFIT
filt_band  = 400;

if f<filt_band %If it is close to zero, just use a lowpass
    [b,a] = butter(3, (f+filt_band/2)*(2/Fs),'low'); %NOTE: if you want to do an 4th order or high filter the tf approach has numerical round off problems and you'll need to use the ss approach comment out below
else
    %This is set up just like the deault filtering in cleandata.m
    half_band = filt_band/2;
    low_freq = max([1,f-half_band]); %lower bound of bandpass. Max is used in case center frequency is less than 1/2 band away from 0 Hz. I use 1 Hz here not 0 Hz, because the filter can't accept 0 Hz
    upp_freq = f+half_band;          %upper bound of bandpass
    %[b,a] = cheby2(3,20,[low_freq upp_freq]*(2/Fs),'bandpass'); %NOTE: if you want to do an 4th order or high filter the tf approach has numerical round off problems and you'll need to use the ss approach comment out below
    [b,a] = butter(3,[low_freq upp_freq]*(2/Fs),'bandpass'); %NOTE: if you want to do an 4th order or high filter the tf approach has numerical round off problems and you'll need to use the ss approach comment out below
end

SDR_filt = filtfilt(b,a,radioin); %zerophase shift
SDR_filt(abs(SDR_filt)<1e-6) = 0; %This was needed in some of our other code to eliminate erroneaouse peaks when no signal is present. Shouldn't be needed here if the signal are know to be there, but leaving it in just in case. 

t_raw = 1/Fs*(0:1:length(radioin)-1)+1/Fs; %create a time vector at raw sample rate. First sample is at 1/Fs
%figure; plot(t_raw,real(radioin));
% hold on; plot(t_raw,real(SDR_filt));

SDR_filt_prod = SDR_filt.*cos(2*pi*f*t_raw');
[b,a] = cheby2(3,20,200*(2/Fs),'low'); %NOTE: if you want to do an 4th order or high filter the tf approach has numerical round off problems and you'll need to use the ss approach comment out below
SDR_shift = filter(b,a,SDR_filt_prod);

% figure; plot(t_raw,abs(SDR_shift));
% hold on; plot(t_raw,real(SDR_filt));

SDR_pulse_sig = abs(SDR_shift);
SDR_pulse_sig = SDR_pulse_sig/max(SDR_pulse_sig);%Normalize

% figure; plot(t_raw,SDR_pulse_sig);
% hold on;plot(t_raw,SDR_pulse_sig>0.10*max(SDR_pulse_sig))

SDR_pulse_logic = SDR_pulse_sig>0.10*max(SDR_pulse_sig);

trans = diff(SDR_pulse_logic);
t_raw(logical(abs(trans)));
% figure; plot(t_raw(1:end-1),trans)

%% GENERATE THE PULSE LOGIC AND TIME DIFFERENCES

%Cleanup oscilating logic
inds_lookback = round(0.01*Fs);
for i = inds_lookback+1:length(trans)
    if all(trans(i-inds_lookback:i-1)==0)
        trans(i) = trans(i);
    else
        trans(i) = 0;
    end
end
%figure; plot(trans);

t_diff_trans = diff(t_raw(logical(abs(trans)))); %Time difference between all transitions

dwell_logic = t_diff_trans>0.5;
durr_logic  = (t_diff_trans>0.01) & (t_diff_trans<0.5);

dwell_list = t_diff_trans(dwell_logic);
durr_list  = t_diff_trans(durr_logic);
smallest_length = min([length(dwell_list),length(durr_list)]);
rep_list = dwell_list(1:smallest_length)+durr_list(1:smallest_length);

durr_mean = mean(durr_list);
rep_mean = mean(rep_list);
durr_std = std(durr_list);
rep_std = std(rep_list);

%% SET OUTPUT VARIABLES
PRI     = rep_mean;
PRI_std = rep_std;
PD      = durr_mean; %Pulse duration
PD_std  = durr_std; %Pulse duration standard deviation
%f;     %carrier frequency difference to from baseband.

end

