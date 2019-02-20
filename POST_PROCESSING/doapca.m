function [out] = doapca(pulse_sig,pulse_yaw,pulse_waypt_num_in,strengthtype,scale,pltcntrl)
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
%                            pulse power         power               power of signal power (power^2) -  this should likely be used. 
%
%   scale               a char array of 'log' or 'linear' to indicate if
%                       the scaling used in the PCA method should be log or
%                       linear
%   pltcntrl            'plot' or 'noplot' to turn on or off the plots
%                       of pulse cloud and bearing estimate
%
%OUTPUTS
%   DOA_calc            a (wx1) numeric vector contain the bearing estimate
%                       from the PCA method from 0-359 degrees with 0 being
%                       the same as the yaw origin (typically N). w here is
%                       the number of waypoints in the pulse_waypt_num 
%                       list, or max(pulse_waypt_num)
%   DOA_tau             a (wx1) numeric vector contain the tau values for
%                       each bearing estimate. w here is the number of
%                       waypoints in the pulse_waypt_num list, or
%                       max(pulse_waypt_num)
%
%
%

num_of_waypts = max(pulse_waypt_num_in);
%If pulse_waypt_num_in is empty, create a vector of all ones so that the
%entire dataset is considered for a single DOA output. If not, use the
%groupings specified in pulse_waypt_num_in

if isempty(pulse_waypt_num_in)
    pulse_waypt_num = ones(size(pulse_sig));
else
    pulse_waypt_num = pulse_waypt_num_in;
end

DOA_calc = zeros(num_of_waypts,1);
DOA_tau = zeros(num_of_waypts,1);


if strcmp(pltcntrl,'ploton')
    figure
    tru_bear_color = [1 0.3 0.3];
    est_bear_color =  [0    0.4470    0.7410];%[0.3 0.5 1];
    mkrsz = 10;
    lnwdth = 1;
end

for i = 1:num_of_waypts
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
        warning('Only %i were detected - insufficient to DOA estimates.', length(curr_pulses))
        wp(2) = 0; wp(1) = 0;tau = 0;
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
    if strcmp(pltcntrl,'ploton')
        subplot(ceil(num_of_waypts/3),min([3,num_of_waypts]),i)
        %%%%%%%%%%polarplot(curr_yaws*pi/180,20*log10(curr_pulses./min(curr_pulses)),'.','Markersize',15,'Color',est_bear_color); hold on;
        %polarplot(curr_yaws*pi/180,(curr_pulses.^2./min(curr_pulses).^2),'.','Markersize',15,'Color',est_bear_color); hold on;
        polarplot(curr_yaws*pi/180,P_all_ang,'.','Markersize',15,'Color',est_bear_color); hold on;
        set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise')

        p_0noise_dir = polarplot([0 atan2(wp(2),wp(1))],line_scale*[0 norm(wp)],'Linewidth',lnwdth,'Color',est_bear_color);hold on
            polarplot([atan2(wp(2),wp(1)) atan2(wp(2),wp(1))-3*pi/180],line_scale*[norm(wp) norm(wp)-norm(wp)*0.1],'b','Linewidth',lnwdth,'Color',est_bear_color)
            polarplot([atan2(wp(2),wp(1)) atan2(wp(2),wp(1))+3*pi/180],line_scale*[norm(wp) norm(wp)-norm(wp)*0.1],'b','Linewidth',lnwdth,'Color',est_bear_color)
        
        max_sig_of_all_wypts = max(pulse_sig(~isnan(pulse_waypt_num))); %look at all the waypoints and find the max pulse amp. exclude pulses not at a waypoint
            max_sig_curr_wypt = max(curr_pulses);
            min_sig_curr_wypt = min(curr_pulses);
        
        %title(['Waypt.',num2str(i),' P_m = ',num2str(20*log10(max_sig_curr_wypt./max_sig_of_all_wypts),2),'dB'],'FontSize',8)
        %title(['\#',num2str(i)],'FontSize',8,'Interpreter','latex')
        
        %%%%%%%%%%set(gca,'RLim',[0 20*log10(max_sig_curr_wypt./min_sig_curr_wypt)])% in dB above min pulse at this waypoint
        %set(gca,'RLim',[0 (max_sig_curr_wypt.^2./min_sig_curr_wypt.^2)])% in dB above min pulse at this waypoint
        set(gca,'Fontsize',12)
        %set(gca,'Fontsize',2)
        %set(gca,'ThetaTickLabel',[])
        set(gca,'Thetatick',0:45:335,'Thetaticklabel',{'N','NE','E','SE','S','SW','W','NW'},'TickLabelInterpreter','tex')
        %set(findall(gcf,'-property','FontSize'),'FontSize',12)
        set(gca,'TickLabelInterpreter','latex')
    end

    DOA_calc(i)  = atan2(wp(2),wp(1));
    DOA_tau(i) = tau;
    %DOA_R(i) = P_all_ang;
end

out{1} = 180/pi*DOA_calc;
out{2} = DOA_tau;

end
