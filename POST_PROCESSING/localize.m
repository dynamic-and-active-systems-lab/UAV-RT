function [position] = localize(X,Y,T,method)
%LOCALIZE estimates a location based on a series of bearing estimates using
%localization methods presented in Lenth's 1972 paper cited below. The 
%methods include MLE with, repeated median regression, and an M-estimator. 
%The outputs is an x-y vector of estimated position found through the
%method requested. 

%NOTE: Lenth uses X (East) and Y (North) position and bearing angles as
%compass headings, which are positive CW from North. This program stays
%consistent with this notation and the X and Y positions, should be the
%East and North positions from which the bearings were measured. This may
%conflict with the frames used by the program calling localize. If the
%X-North Y-East coordinate frame is used, x and y positions entered to this
%function can simply be switched, and the output positions should also then
%be switched. 

% INPUTS:
%X           -    n x 1 vector of x-location where pulse_sig was measured. 
%                       X here is distance East "
%Y           -    n x 1 vector of y-location where pulse_sig was measured
%T           -    n x 1 vector of bearing angles in DEGREES from North with
%                       positive in CW direction. Standard compass bearing.
%method      -    a char array of selected method:
%                   Options are 'MLE', 'RMR', or 'MEST'

% OUTPUTS:
%position    -    1 x 2 estimate of [x,y] location of signal origin


%Methods implemented here are from 
%Lenth, Russell. On Finding the Source of a Signal. Technometeric. Vol. 23,
%No. 2, 1981. 

%MLE: Maximum likely hood esitmator
%RMR: Repeated median regression
%MEST: M-estimation with Andrews psi function with c = 1.5

%Make sure these are all column vectors
if ~iscolumn(X)
    X = X';
end
if ~iscolumn(Y)
    Y = Y';
end
if ~iscolumn(T)
    T= T';
end


max_iter = 1000;
if ~isequal(numel(X),numel(Y),numel(T)) %Check to see that vectors are same length
    error('Number of elements of inputs must be equal')
else
    %Eliminate any that have NaN X, Y, or T values. 
    nan_mask = ~isnan(X)&~isnan(Y)&~isnan(T);
    X = X(nan_mask);
    Y = Y(nan_mask);
    T = T(nan_mask);
    num_bears = length(X);
end

%Standardize the shape of all the input vectors
X = reshape(X,num_bears,1);
Y = reshape(Y,num_bears,1);
T = reshape(T,num_bears,1);

%Prep some variables that will be needed for RMR
if strcmp(method,'RMR') 
    %Create the element selection vectors
    %This creates the list of all combintations
    bc = nchoosek(1:num_bears,2);
    %Create the element selection vectors
    %This creates the list of all combintations
    bc = nchoosek(1:num_bears,2);
    bcA = bc(:,1);
    bcB = bc(:,2);
    %Generate the positions and angles of the different selected elements for
    %comparison: A vs B.
    XA = X(bcA);
    XB = X(bcB);
    YA = Y(bcA);
    YB = Y(bcB);
    TA = T(bcA);
    TB = T(bcB);
    %Find the intersection points for all the A and B lines
    X_int = (YA - YB + XB.*cotd(TB) - XA.*cotd(TA))./(cotd(TB)-cotd(TA));
    Y_int = cotd(TA).*X_int+YA-XA.*cotd(TA);
    %Some test cases had a few NaN locations, so only track those that have no
    %NaNs either X_int and Y_int. We redefine everything here to rid the record
    %of those points.
    valid_pos_msk = ~isnan(X_int)&~isnan(Y_int);
    X_int = X_int(valid_pos_msk);
    Y_int = Y_int(valid_pos_msk);  
end



%% MAXIMUM LIKELYHOOD ESTIMATOR

theta_i = pi/180*(90-T);%On page 1, paragraph of the Lenth paper. 
xi = X;% Lenth says x is E and y is N, but we use x N and y E, 
yi = Y;
si = sin(theta_i);
ci = cos(theta_i);
zi = si.*xi-ci.*yi;
pos(:,1) = [sum(si.*si) -sum(ci.*si);-sum(si.*ci) sum(ci.*ci)]\[sum(si.*zi);-sum(ci.*zi)];

if  strcmp(method,'MLE')
    i = 1;res = 1;
    while res(i)>0.001 && i<max_iter
        x = pos(1,i);
        y = pos(2,i);
        di = sqrt((x-xi).^2+(y-yi).^2);
        si_str = (y-yi)./di.^3;
        ci_str = (x-xi)./di.^3;
        
        pos(:,i+1) = [sum(si.*si_str) -sum(ci.*si_str);-sum(si.*ci_str) sum(ci.*ci_str)]\[sum(si_str.*zi);-sum(ci_str.*zi)];
        i = i+1;
        res(i) = norm(pos(:,i)-pos(:,i-1));
    end
    if i==max_iter
        warning(['MLE method did not converge after ',num2str(max_iter),' iterations.'])
        pos_MLE = [NaN;NaN];
    else
        pos_MLE = pos(:,end);
    end

    %disp(['loc_error_MLE = ',num2str(norm(pos_MLE))])
    position = pos_MLE;
    
%% REPEATED MEDIAN REGRESSION



elseif strcmp(method,'RMR')
    
    %
    for i = 1:num_bears
        tick = 1;
        for j = [1:i-1,i+1:num_bears]
            ind_sel = find( (bc(:,1) == i  & bc(:,2) == j) | (bc(:,2) == i  & bc(:,1) == j)); %Find the entry that had both i and j in either the first or second columns
            int_log(:,tick) = [X_int(ind_sel);Y_int(ind_sel)]; %create a log of all i-th ray intersecting with all other rays
            tick = tick+1;
        end
        MED_INT_LOG(:,i) = median(int_log,2);%take the mean of the i-th ray's intersection locations.
    end
    
    pos_RMR  = median(MED_INT_LOG,2);
    %disp(['loc_error_RMR = ',num2str(norm(pos_RMR))])
    position = pos_RMR;
    

%% M-ESTIMATOR
elseif strcmp(method,'MEST')   
    theta_i = pi/180*(90-T);
    xi = X;
    yi = Y;
    si = sin(theta_i);
    ci = cos(theta_i);
    zi = si.*xi-ci.*yi;
    
    si_hat_str = si;
    ci_hat_str = ci;
    
    wi = ones(size(X));
    pos3(:,1) = ([sum(wi.*si.*si_hat_str) -sum(wi.*ci.*si_hat_str);-sum(wi.*si.*ci_hat_str) sum(wi.*ci.*ci_hat_str)])\[sum(wi.*si_hat_str.*zi);-sum(wi.*ci_hat_str.*zi)];
    
    i = 1;clear res;res(1) = 1;
    while res(i)>0.001 && i<max_iter
        x = pos3(1,i);
        y = pos3(2,i);
        
        mu_i = atan2(y-yi,x-xi);
        c = 1.5;
        phi = theta_i-mu_i;
        
        C_bar_w = sum(wi.*cos(theta_i-mu_i))/sum(wi);
        kappa = 1/(2*(1-C_bar_w)+(1-C_bar_w)^2*(0.48794-0.82905*C_bar_w-1.3915*C_bar_w^2)/C_bar_w);
        %Without kappa approximation
        %     bess_input = [C_bar_w-0.0001,C_bar_w];%Create input to bessel function to take derivative; Small step back that we'll use to approx derivative.
        %     A = diff(log(besseli(0,bess_input)))/diff(bess_input);
        %     kappa = A^-1;

        t = sqrt(2*kappa*(1-cos(phi)));
        psi = c*sin(t/c).*(abs(t)<c*pi);
        wi = psi./t;
        
        di = sqrt((x-xi).^2+(y-yi).^2);
        si_hat_str = (y-yi)./di.^3;
        ci_hat_str = (x-xi)./di.^3;
        
        pos3(:,i+1) = ([sum(wi.*si.*si_hat_str) -sum(wi.*ci.*si_hat_str);-sum(wi.*si.*ci_hat_str) sum(wi.*ci.*ci_hat_str)])\[sum(wi.*si_hat_str.*zi);-sum(wi.*ci_hat_str.*zi)];
        %([sum(si.*si_str) -sum(ci.*si_str);-sum(si.*ci_str) sum(ci.*ci_str)])\[sum(si_str.*zi);-sum(ci_str.*zi)];
        i = i+1;
        res(i) = norm(pos3(:,i)-pos3(:,i-1));
        
    end
	if i==max_iter
        warning(['MLE method did not converge after ',num2str(max_iter),' iterations.'])
        pos_MEST = [NaN;NaN];
    else
        pos_MEST = pos3(:,end);
    end
    %disp(['loc_error_MEST = ',num2str(norm(pos_MEST))])
    position = pos_MEST;
    
else
    warning('Method entered is not supported')
    position = [NaN, NaN];
end

end

