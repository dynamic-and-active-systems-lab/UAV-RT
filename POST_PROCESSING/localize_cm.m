function [position,intmap] = localize_cm(X,Y,T,pulse_sig,tau,tau_lim,ang_sep_lim)
%LOCALIZE_CM estimates a location based on a series of bearing estimates 
%using the center of mass localization methods. In this method the
%intersection of all bearings are found and then weighed by the product of
%the signal strengths that created them. The location is estimated by
%finding the center of mass of these weightings
%

%NOTE: In order to stay consistent with LOCALIZE.M this function uses the
%same reference frame used by Lenth (1972) cited below. 
%Lenth uses X (East) and Y (North) position and bearing angles as
%compass headings, which are positive CW from North. This program stays
%consistent with this notation and the X and Y positions, should be the
%East and North positions from which the bearings were measured. This may
%conflict with the frames used by the program calling localize. If the
%X-North Y-East coordinate frame is used, x and y positions entered to this
%function can simply be switched, and the output positions should also then
%be switched. If the ang_sep_lim or the tau_lim limits are not met, the
%function still provides a localization esitmate, but warns the user that
%there may be significant localization error. 

% INPUTS:
%X           -    n x 1 vector of x-location where pulse_sig was measured. 
%                       X here is distance East "
%Y           -    n x 1 vector of y-location where pulse_sig was measured
%T           -    n x 1 vector of bearing angles in DEGREES from North with
%                       positive in CW direction. Standard compass bearing.
%pulse_sig   -    n x 1 vector of signal strength for each of the bearing
%                       estimates. Could be average or maximum power or
%                       signal amplitude.
%tau         -    n x 1 vector of tau values of the DOA estimates
%tau_lim     -    1 x 1 scalar of the minimum tau to consider
%ang_sep_lim -    1 x 1 scalar of the minimum angular separation between 
%                       bearing to consider valid when developing 
%                       intersection points in degrees. 

%OUTPUTS:
%position    -    1 x 2 estimate of [x,y] location of signal origin
%intmap      -    k x 5 matrix of x,y position (columns 1 and 2) of each 
%                       bearing intersection, the weighting of that point
%                       (column 3) and the bearing numbers that
%                       created the intersection (columns 4 and5)

%Coordinate system implemented here is consistent with that from 
%Lenth, Russell. On Finding the Source of a Signal. Technometeric. Vol. 23,
%No. 2, 1981. 

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
if ~iscolumn(pulse_sig)
    pulse_sig= pulse_sig';
end
if ~iscolumn(tau)
    tau= tau';
end




if ~isequal(numel(X),numel(Y),numel(T),numel(pulse_sig),numel(tau)) %Check to see that vectors are same length
    error('Number of elements of inputs must be equal')
else
   %Eliminate any that have NaN X, Y, or T values. 
    nan_mask = ~isnan(X)&~isnan(Y)&~isnan(T);
    X = X(nan_mask);
    Y = Y(nan_mask);
    T = T(nan_mask);
    pulse_sig = pulse_sig(nan_mask);
    tau = tau(nan_mask);
    num_bears = length(X);
end

%Standardize the shape of all the input vectors
pulse_sig = reshape(pulse_sig,num_bears,1);
X = reshape(X,num_bears,1);
Y = reshape(Y,num_bears,1);
T = reshape(T,num_bears,1);
tau = reshape(tau,num_bears,1);

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
XA = XA(valid_pos_msk);
XB = XB(valid_pos_msk);
YA = YA(valid_pos_msk);
YB = YB(valid_pos_msk);
TA = TA(valid_pos_msk);
TB = TB(valid_pos_msk);
bcA = bcA(valid_pos_msk);
bcB = bcB(valid_pos_msk);



%Here is the product of the pulse signals. Use this to weight the
%localization
pulse_prod = pulse_sig(bcA).*pulse_sig(bcB);


%R = sqrt(X_int.^2+Y_int.^2);

%Calculate the lists of bearing angle to the interection points from the
%locations of the parent line origins. The develop a clip element when the
%difference between these bearings is less than the threshold. 
T_INT_A = 180/pi*atan2(X_int-XA,Y_int-YA);
T_INT_B = 180/pi*atan2(X_int-XB,Y_int-YB);
sep_list = (abs(wrapTo360(T_INT_B)-wrapTo360(T_INT_A)));
clip_ang_sep = sep_list>ang_sep_lim;
if all(~clip_ang_sep) %if all clip are zero
    warndlg(['Bearing estimate threshold not met. ' ,...
              'No location estimate generated. ',...
              'Try separating waypoints, or decreasing threshold. ',...
              'The [min max] angular separation of your bearings is [',...
              num2str(min(sep_list)),', ',...
              num2str(max(sep_list)),'] degrees. REPORTED LOCATION ESTIMATE ',...
              'MAY HAVE SIGNIFICANT ERROR.']);
	%Still give provide an estimate, because there is no other option
	clip_ang_sep = ones(size(clip_ang_sep));
end

%This is a simple check to see if the bearing to the intersection point is
%in front of or behind the waypoint. If it is in front, the theta value for
%that waypoint should be the same as the intersection point. If the
%intersection occurs behind the waypoint, the angle to the intersection
%point will be 180+T_A, so TA and T_INT_A will not be equal. We use the
%round command because roundin error in the generation of T_INT_A will
%result in small differences between TA and T_INT_A.
frontA = round(TA)==round(T_INT_A);
frontB = round(TB)==round(T_INT_B);
front = frontA|frontB;
clip_front = front;


%clip = clip_front;
%clip_range = (sqrt((mean(X)-X_int).^2+(mean(Y)-Y_int).^2))<1000;%Range from the average vehilce position to all intersections.  R<100;
%clip_pulse_prod = pulse_prod>0.5*max(pulse_prod);
%clip = clip_front&clip_range&clip_pulse_prod;

%Create a clip element that only allows bearing where the TAU value was
%above the threshold. Need to filter out both the A and B lines that might
%have had tau values below the threshold.
clip_tau = (tau(bcA)>=tau_lim)&(tau(bcB)>=tau_lim);
if sum(clip_tau)<2 && length(bcA)>1%if fewer than 2 good taus and we had more than one combination of points...
    warndlg(['Did not have sufficient tau values for least two bearings. ' ,...
              'May need better data or try decreasing tau threshold. ',...
              'The [min max] tau of your bearings is [',...
              num2str(min(tau)),', ',...
              num2str(max(tau)),']. REPORTED LOCATION ESTIMATE ',...
              'USES ALL BEARINGS AND MAY HAVE SIGNIFICANT ERROR.']);
	%Still give provide an estimate, because there is no other option
    clip_tau = ones(size(clip_tau));
end

%Combine all the clip elements
clip = clip_front&clip_tau&clip_ang_sep;


%Generate the be belief vector (zero out the ones we chose to clip)
belief = pulse_prod.*clip;

% figure;
% plot3(X_int(clip),Y_int(clip),pulse_prod(clip)/max(pulse_prod(clip))*500,'.')
% hold on
% plot3(XA,YA,zeros(size(XA)),'r.','Markersize',10)
% 
% figure;
% plot(X_int(clip),Y_int(clip),'.')
% hold on
% plot(X,Y,'r.','Markersize',10)
% axis equal; grid on


%Use the CM method to esitmate the location. 
X_loc = 1/sum(pulse_prod(clip))*sum(pulse_prod(clip).*X_int(clip));
Y_loc = 1/sum(pulse_prod(clip))*sum(pulse_prod(clip).*Y_int(clip));

%loc_error = sqrt(X_loc^2+Y_loc^2)

% [xi,yi] = meshgrid(-1000:1:10000, -1000:1:1000);
% zi = griddata(X_int(clip),Y_int(clip),pulse_prod(clip),xi,yi);
% figure
% surf(xi,yi,zi,'Edgecolor','none');
% hold on;
% stem3(X_loc,Y_loc,2e-6,'r')

pos_CM = [X_loc;Y_loc];
%disp(['loc_error_CM = ',num2str(norm(pos_CM))])


position = pos_CM;
intmap = [X_int, Y_int, belief, bcA, bcB];


end

