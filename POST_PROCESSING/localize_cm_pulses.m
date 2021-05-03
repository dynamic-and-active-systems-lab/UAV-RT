function [position,intmap] = localize_cm_pulses(X,Y,pulse_sig)
%LOCALIZE_CM estimates a location based on a series of pulses
%using the center of mass localization methods. In thismethod each pulse is
%then weighed by itssignal strengths. The location is estimated by
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
%be switched. 

% INPUTS:
%X           -    n x 1 vector of x-location where pulse_sig was measured. 
%                       X here is distance East "
%Y           -    n x 1 vector of y-location where pulse_sig was measured
%pulse_sig   -    n x 1 vector of signal strength for each of the bearing
%                       estimates. Could be average or maximum power or
%                       signal amplitude.

%OUTPUTS:
%position    -    1 x 2 estimate of [x,y] location of signal origin

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
if ~iscolumn(pulse_sig)
    pulse_sig= pulse_sig';
end



if ~isequal(numel(X),numel(Y),numel(pulse_sig)) %Check to see that vectors are same length
    error('Number of elements of inputs must be equal')
else
   %Eliminate any that have NaN X, Y, or T values. 
    nan_mask = ~isnan(X)&~isnan(Y);
    X = X(nan_mask);
    Y = Y(nan_mask);
    pulse_sig = pulse_sig(nan_mask);
    num_pulses = length(X);
end

%Standardize the shape of all the input vectors
pulse_sig = reshape(pulse_sig,num_pulses,1);
X = reshape(X,num_pulses,1);
Y = reshape(Y,num_pulses,1);


%Find the points for all points
%These are redefine just to stay consisten with the original code this was
%based on (localize_cm.m)
X_int = X;
Y_int = Y;

%Some test cases had a few NaN locations, so only track those that have no
%NaNs either X_int and Y_int. We redefine everything here to rid the record
%of those points. 
valid_pos_msk = ~isnan(X_int)&~isnan(Y_int);
X_int = X_int(valid_pos_msk);
Y_int = Y_int(valid_pos_msk);



%Use the CM method to esitmate the location. 
X_loc = 1/sum(pulse_sig)*sum(pulse_sig.*X_int);
Y_loc = 1/sum(pulse_sig)*sum(pulse_sig.*Y_int);

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
intmap = [X_int, Y_int, pulse_sig];


end

