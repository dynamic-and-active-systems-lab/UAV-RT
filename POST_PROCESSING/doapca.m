function [out] = doapca(pulse_sig,pulse_yaw,pulse_waypt_num_in,total_waypts,strengthtype,scale,pltcontainer)
%DOAPCA developes a bearing estimate for a series of received radio pulses
%based on the principle component analysis method. 
%   This function conductes a principle component analysis type bearin
%   estimated on signal pulses in the pulse_sig vector. Each pulse signal
%   corresponds to a given heading (pulse_yaw) and a waypoint number
%   (pulse_waypt_num). If pulse_waypt_num is empty, the analysis is run on
%   the entire dataset. If pulse_waypt_num indicates that some pulses are
%   grouped at certain waypoints, the function will run the PCA on the
%   individual waypoints and report bearing estimates for each. 
%
%
%
%INPUTS
%   pulse_sig           a (px1) numeric vector with signal amplitude for 
%                       each pulse. p is the number of pulses
%   pulse_yaw           a (px1) numeric vector with the heading when the
%                       pulse was received in degrees
%   pulse_waypt_num_in  a (px1) numeric vector with the waypoint number of
%                       each pulse and NaNs for pulses not associated with
%                       a waypoint. Enter [] for this variable if you want
%                       to consider the entire dataset. 
%   total_waypts        an integer that specifies the total number of
%                       waypoints of this flight. This is used to ensure
%                       the output contains the same number of elements as
%                       the number of waypoints. max(pulse_waypt_num_in)
%                       may be less than total_waypoints, if some waypoints
%                       did not contain pulses. Outputs for waypoints that
%                       don't have pulse data contain NaN
%   strengthtype        a char array of 'power' or 'amplitude' indicating
%                       the strength used in the PCA method should use the
%                       signal amplitude or power. If power is used, the
%                       square of the inputs pulse signal is used. If
%                       amplitude if fed as pulse_sig and power is
%                       selected, the inputs to the PCA is power. If power
%                       is fed as pulse_sig and amplitude is selected as
%                       strength type, the input to the PCA is power. See
%                       table below
%                            pulse_sig           strengthtype        PCA input
%                       ------------------------------------------------------
%                            pulse amplitudes    amplitude           signal amplitude
%                            pulse amplitudes    power               signal power
%                            pulse power         amplitude           signal power
%                            pulse power         power               power of signal power (power^2) -  this shouldn't likely be used. 
%
%   scale               a char array of 'log' or 'linear' to indicate if
%                       the scaling used in the PCA method should be log or
%                       linear
%   pltcontainer        This is a handle of a container (figure handle,
%                       etc.) where the bearing estimate figures should go.
%                       This function generates a set of subplots within
%                       this this container, one for each waypoint that
%                       shows pulse cloud and bearing estimate of each
%                       waypoint. If you don't want anything plotted, use
%                       the empty values [] here. 
%
%OUTPUTS
%   DOA_calc            a (wx1) numeric vector contain the bearing estimate
%                       from the PCA method from 0-359 degrees with 0 being
%                       the same as the yaw origin (typically N). w here is
%                       the number of waypoints in the pulse_waypt_num 
%                       list, or max(pulse_waypt_num). If waypoints were
%                       skipped, NaN is reported for DOA_calc. 
%   DOA_tau             a (wx1) numeric vector contain the tau values for
%                       each bearing estimate. w here is the number of
%                       waypoints in the pulse_waypt_num list, or
%                       max(pulse_waypt_num).If waypoints were
%                       skipped, NaN is reported for DOA_tau. 
%
%
%

num_of_waypts = total_waypts; max(pulse_waypt_num_in);

%If pulse_waypt_num_in is empty, create a vector of all ones so that the
%entire dataset is considered for a single DOA output. If not, use the
%groupings specified in pulse_waypt_num_in

if isempty(pulse_waypt_num_in)
    pulse_waypt_num = ones(size(pulse_sig));
else
    pulse_waypt_num = pulse_waypt_num_in;
end

%DOA_calc = zeros(num_of_waypts,1);
%DOA_tau = zeros(num_of_waypts,1);
%Pre allocate so that if we skip over a waypoint, the output will be NaN
%for that waypoint
DOA_calc = NaN(num_of_waypts,1);
DOA_tau = NaN(num_of_waypts,1);

%This is the list of waypoints that should be evaluted for DOA.
%Eliminate the NaN values, and create a unique list
waypts_with_pulses = unique(pulse_waypt_num_in(~isnan(pulse_waypt_num_in)));

if ~isempty(pltcontainer)
    %figure
    %tru_bear_color = [1 0.3 0.3];
    est_bear_color =  [0    0.4470    0.7410];%[0.3 0.5 1];
    mkrsz = 10;
    lnwdth = 1;
end

for i = waypts_with_pulses%1:num_of_waypts
    curr_pulses = pulse_sig(pulse_waypt_num == i);
    curr_yaws = pulse_yaw(pulse_waypt_num == i);
    
    if strcmp(strengthtype,'power')
        P_all_ang_unscaled = (curr_pulses.^2./min(curr_pulses).^2)';
    elseif strcmp(strengthtype,'amplitude')
        P_all_ang_unscaled = (curr_pulses./min(curr_pulses))';
    end
    
    if strcmp(scale,'linear')
        P_all_ang = P_all_ang_unscaled;
    elseif strcmp(scale,'log')
        P_all_ang = log10(P_all_ang_unscaled);
    end
    
    angs = curr_yaws*pi/180;
    
    if length(curr_pulses)<4
        warning(['Only ',num2str(length(curr_pulses)),' pulse(s) detected at waypoint #',num2str(i),' - insufficient to make DOA estimate.'])
        wp(2) = NaN; wp(1) = NaN;tau = NaN; line_scale = 0;
    else
	    Pe_star_dB = [P_all_ang.*cos(angs),P_all_ang.*sin(angs)];
        n = length(Pe_star_dB);
        Pe_dB = (eye(n)-1/n*ones(n))*Pe_star_dB;
        Pavg = 1/n*Pe_star_dB'*ones(n,1);
        [~, SdB, VdB] = svd(Pe_dB);
        w1 = VdB(:,1);
        w2 = VdB(:,2);
        wp = Pavg'*w1/(norm(Pavg)*norm(w1))*w1;
        beta = norm(Pavg)^2/SdB(1,1)^2;
        tau = 1-SdB(2,2)^2/SdB(1,1)^2;
        line_scale = max(P_all_ang)/norm(wp);%the wp size changes if w1
       
    end
    if ~isempty(pltcontainer)
        axcurr = subplot(ceil(num_of_waypts/3),min([3,num_of_waypts]),i,polaraxes,'Parent',pltcontainer);        
        polarplot(axcurr,curr_yaws*pi/180,P_all_ang,'.','Markersize',15,'Color',est_bear_color); 
        hold(axcurr,'on');
        set(axcurr,'ThetaZeroLocation','top','ThetaDir','clockwise')

        p_0noise_dir = polarplot(axcurr,[0 atan2(wp(2),wp(1))],line_scale*[0 norm(wp)],'Linewidth',lnwdth,'Color',est_bear_color);hold on
            polarplot(axcurr,[atan2(wp(2),wp(1)) atan2(wp(2),wp(1))-3*pi/180],line_scale*[norm(wp) norm(wp)-norm(wp)*0.1],'b','Linewidth',lnwdth,'Color',est_bear_color)
            polarplot(axcurr,[atan2(wp(2),wp(1)) atan2(wp(2),wp(1))+3*pi/180],line_scale*[norm(wp) norm(wp)-norm(wp)*0.1],'b','Linewidth',lnwdth,'Color',est_bear_color)
        
        max_sig_of_all_wypts = max(pulse_sig(~isnan(pulse_waypt_num))); %look at all the waypoints and find the max pulse amp. exclude pulses not at a waypoint
            max_sig_curr_wypt = max(curr_pulses);
            min_sig_curr_wypt = min(curr_pulses);
        
        set(axcurr,'Fontsize',12)
        set(axcurr,'Thetatick',0:45:335,'Thetaticklabel',{'N','NE','E','SE','S','SW','W','NW'},'TickLabelInterpreter','tex')
%        set(gca,'TickLabelInterpreter','latex')
    end

    DOA_calc(i)  = atan2(wp(2),wp(1));
    DOA_tau(i) = tau;
end

out{1} = 180/pi*DOA_calc;
out{2} = DOA_tau;

end
