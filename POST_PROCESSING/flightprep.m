function [a_sum, waypt, home] = flightprep(data_filename,varargin)
%FLIGHT_PREP reads in a flight log and preps is for use by other functions.
%   This function will read in a file and create the flight log matrix
%   necessary by other functions. I will convert the lat long to x and y
%   distances from the home location. The default home location is first
%   lat lon in the flight log. Options 'lat' and 'lon' can be used to
%   specift an alternate home location for the (0,0) position. 
%   
%   Flight log data must have first column for time, and latitude and 
%   longitude in other columns. Other columns are accepted, but not 
%   required.
%   The output frame for x and y is North and East, respectively so that
%   positive heading angle and x/y positions are consistent. 
%
%Inputs:
%   data_filename       string filename or path to flight log
%                       data. Must be a tab delimited file with column headings
%                         Time_s
%                         Roll_rad; %roll in radians
%                         Pitch_rad; %pitch in radians
%                         Yaw_rad; %yaw in radians
%                         Alt_m_AGL; %alt in m AGL
%                         Lat; %Latitude in decimal forms
%                         Long; %Longitude in decimal form in range[-180 180]
%   
%Options (don't need to spell out the name if value comes right after data_filename)
%   'timeselect'       [t_start t_end] start time and end time that should
%                       be considered, measured from the start of the 
%                       record. Enter 'all' for entire time.
%   'plot'              Turns on plotting results. Default is 'noplot'. Use
%                       'plot' to plot vehicle state results
%
%Parameters (need to spell out the name, then enter the value)
%   'latlon'             [latitude longitude] of home location 
%                       (Launch location used if not entered
%   'fs'                sampling rate of telemetry in Hz. Default is 10 if
%                       not entered
%   'waypt_vel'         velocity (m/s) threshold below which you can be 
%                       considered at a waypoint. Default is 0.5 m/s if
%                       not entered.
%   'waypt_time'        min and max time required to be below waypt_vel to 
%                       define a waypoint. Default is [0 10000] s if not entered.
%   'waypt_alt'         This is the waypoint altitudes (AGL) definition
%                       listting. If entered as a 1x1 the waypoints is 
%                       considered withing +/-10% of that altitudes. 
%                       If not entered +/-10% if max alt is used. 
%                       If listed as an [nx2], each row vector represents
%                       the bounds of a waypoint definition. For example if
%                       one waypoint in the flight was somewhere between 
%                       20 and 22 m, and another was between 45 and 46 m,
%                       the entry would be [20 22; 45 46]. In this format,
%                       the lower altitude must be in the first column. 
%   'axis_handles'      Handle of axis where data should be plotted. If
%                       empty, function will create a new figure. 



%Outputs:
%   a_sum               numeric array of the flight data record 
%                       Definitions of the a_sum matrix columns
%                           1  delta t in sec from start
%                           2  roll (deg)
%                           3  pitch(deg),
%                           4  total yaw (deg from true N, CW+), 
%                                Note that this is the unwrapped angle.
%                                Wrapping to 0-360 can cause problems with
%                                interpolations. If the vehilce starts at
%                                45 deg and spins 5 times, this column
%                                would start at 45 and end at 45+ 360*5 =
%                                1845 deg. 
%                           5  alt (m), 
%                           6  latitude, 
%                           7  longitude, 
%                           8  x pos (m) with positive +N
%                           9  y pos (m) with positive +E
%                           10 at waypoint logic
%                           11 at altitude logic
%   waypt               numeric array of waypoint summary information
%                       Definitions of the waypt matrix columns
%                           1  mean x position (m) of waypt, 
%                           2  mean y position (m) of waypt
%                           3  mean altx position (m) of waypt,
%                           4  time at start of waypoint 
%                           5  time at end of waypoint
%   home                The [latitude, longitude] of the home position, as
%                       specified by the use or if not entered, it is taken
%                       to be the lat lon at the beginning of the data
%                       record
%


%Examples
% flight matrix = flight_prep('flight_data.txt')
%   Read the flight_data.txt log and creates outputs using default options
%   and parameters for the entire dataset
%
% flight matrix = flight_prep('flight_data.txt',[5 50],'lat',37.12345,'lon',-117.12345)
%   Read the flight_data.txt log and creates outputs using the specified 
%   lat lon as the (0,0) position and only considers time from 5 to 50 s. 



%% Parse the input variables
%Modeled from 
% https://www.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
p = inputParser;

default_latlon = [NaN NaN];%Set to Nan. We'll catch this later and enter default values from the data record after it is loaded.
default_fs = 10;
default_v_thresh_at_wypt = 0.5; 
default_t_bounds = 'all'; %Set to zero so we can use the bounds of the time vector later. 
default_waypt_dwell = 'auto';
default_waypt_alt = NaN;%Set to Nan. We'll catch this later and enter default values from the data record after it is loaded.
default_plot = 'noplot';%Don't plot by default
default_axis = [];   %Create empty variable for figure handle. 
    

%This checks to see if the time entered is either 'all' or an array of 1x2
%and that time is increasing. 
checktime = @(x) strcmp(x,'all')||(isnumeric(x) && size(x,1)==1 && size(x,2)==2 && x(2)>x(1));
checkdwell = @(x) strcmp(x,'auto')||(isnumeric(x) && size(x,1)==1 && size(x,2)==2 && x(2)>x(1));
allhandle = @(x) all(ishandle(x))||isempty(x); %Check all array elements are handles or is empty for default case

addRequired(p,'file_location',@ischar);
addOptional(p,'timeselect',default_t_bounds,checktime);
addOptional(p,'plot',default_plot,@ischar);
addParameter(p,'latlon',default_latlon,@isnumeric);
addParameter(p,'fs',default_fs,@isnumeric);
addParameter(p,'waypt_vel',default_v_thresh_at_wypt,@isnumeric);
addParameter(p,'waypt_time',default_waypt_dwell,checkdwell);
addParameter(p,'waypt_alt_defs',default_waypt_alt,@isnumeric);
addParameter(p,'axis_handles',default_axis,allhandle);%If array, all must be a handle


parse(p,data_filename,varargin{:})

telemfile = p.Results.file_location;
lat_home = p.Results.latlon(1);
lon_home = p.Results.latlon(1);
fs_telem = p.Results.fs;
v_thresh_at_wypt = p.Results.waypt_vel;
t_bounds = p.Results.timeselect; %This could be a string 'all' or the bounds of time [t_start, t_end]
waypt_dwell = p.Results.waypt_time;
waypt_alt = p.Results.waypt_alt_defs;
plot_control =  p.Results.plot;
ax2plot =  p.Results.axis_handles;

%% Get plot axes ready
if strcmp(plot_control,'plot')
    if isempty(ax2plot)
        figure;
        subplot(2,1,1) %plot the x-y trac
        ax2plot(1) = gca;
        subplot(2,1,2) %plot the x-y trac
        ax2plot(2) = gca;
    end
end       

%% Read Flight Data
A_telem = tdfread(telemfile);
a_sum(:,1)=A_telem.Time_s-A_telem.Time_s(1);
a_sum(:,2)=rad2deg(A_telem.Roll_rad); %roll
a_sum(:,3)=rad2deg(A_telem.Pitch_rad); %pitch
%a_sum(:,4)=wrapTo360(rad2deg(A_telem.Yaw_rad)); %yaw correcton to go from 0 to 360 vs -180 to 180
a_sum(:,4)=180/pi*unwrap(wrapTo2Pi(A_telem.Yaw_rad)); %yaw in degrees from start. We unwrap, because we have to lower the resolution below to clean up repeated time logs, and if we wrap to 2pi here, we would end up with yaw records about halfway between 2pi and 0. We'll rewrap to 360 after the code below.
a_sum(:,5)=A_telem.Alt_m_AGL; %alt
a_sum(:,6)=A_telem.Lat; %Lat
a_sum(:,7)=A_telem.Long; %Long


%%
%This section of code cleans up the position data that seems to report the
%same position with subsequent time stamps, may because of the position
%resolution. We first lower the time resolution by about a factor of 2 and
%then resample at the original time scale. The result though creates a Nan
%as the last data point because of the interpolation. To deal with this, we
%just delete the last data point.
time        = a_sum(:,1);
time_lowres = time(1):0.2:time(end);

a_sum_lowres = interp1(time,a_sum,time_lowres);
a_sum = interp1(time_lowres,a_sum_lowres,time); %redefine a_sum 

%a_sum(:,4) = wrapTo360(180/pi*a_sum(:,4)); %No longer necessary as we are using the total rotation angle in a_sum(:,4) rather than the wrapped heading angle. Was...Convert to degrees and wrap to 360
%Clip the last second of points. 
n_1_sec = fs_telem;
a_sum = a_sum(1:end-n_1_sec,:);
time  = time(1:end-n_1_sec,:);


%%
% Now that entries aren't repeated, lets calculate x and y positions
%relative to take off position calculator (x and y)
lats=deg2rad(a_sum(:,6));
lons=deg2rad(a_sum(:,7));

%If not specified use the first location as home.
%lat_home and lon_home are already defined in the parse if specified by 
%the user so no else statements is needed here. 
if isnan(lat_home)||isnan(lon_home)
    lat_home =lats(1);
    lon_home =lons(1); 
end

dlons = lons-lon_home;
dlats = lats-lat_home;
b = (sin(dlats/2)).^2 + cos(lat_home) .* cos(lats) .* (sin(dlons/2)).^2;
c = 2 * atan2(sqrt(b),sqrt(1-b));
R = 6371000 ;%R is the radius of the Earth
d = R .* c; %diameter
B = atan2(cos(lat_home).*sin(lats)-sin(lat_home)*cos(lats).*cos(lons-lon_home),sin(lons-lon_home).*cos(lats)); %angle
% %USING A FRAME WITH +X = EAST AND +Y = NORTH
% a_sum(:,8)=d.*cos(B); %This is distance in the East Direction in meters
% a_sum(:,9)=d.*sin(B); %This is distance in the North direction in meters
%USING A FRAME WITH +X = NORTH AND +Y = EAST
a_sum(:,8)=d.*sin(B); %This is distance in the North direction in meters = X Position
a_sum(:,9)=d.*cos(B); %This is distance in the East Direction in meters = Y Position



%Now set the waypt alt if provided by the user 
if isnan(waypt_alt)
    waypt_alt = max(A_telem.Alt_m_AGL)/1.1;%We use +/- 10% of this values later so range of consideration would be [0.9*max/1.1 =  0.82*max, 1.1*max/1.1 = 1*max]
end

%Create some variables to make code more readable. 
x_pos       = a_sum(:,8);%This is distance in the North direction in meters
y_pos       = a_sum(:,9);%This is distance in the East direction in meters
alt         = a_sum(:,5);
yaw_ang     = a_sum(:,4);


%This block of code calculates a moving mean of the velocity to determine
%when we are at a way point. 
    %R = sqrt(x_pos.^2+y_pos.^2);
    %R_dot = diff(R)./diff(time);
    x_dot = diff(x_pos)./diff(time);
    y_dot = diff(y_pos)./diff(time);
    V = sqrt(x_dot.^2+y_dot.^2);
    %Number of points for 3 sec
    n_3_sec = 3*n_1_sec;
    V_smooth = movmean(V,n_3_sec);%dt is 0.1s, so k = 10 is 1 s

%Limit time search to those provided by user or default (entire time
%record)
if strcmp(t_bounds,'all') %if t_bounds is a vector, this will return logical 0 and use entire dataset
    t_start = time(1);
    t_end = time(end);
else
    t_start = max([t_bounds(1),time(1)]); %just in case entered value is less than the start time.
    t_end = min([t_bounds(2),time(end)]);%just in case enter value is higher than the start time. 
end

%Create the logic for when we are at the prescribed altitudes.
msk_at_alt = zeros(size(time));
if numel(waypt_alt)==1 %Single altitude given
    msk_at_alt = (alt>=0.9*waypt_alt&...
                  alt<=1.1*waypt_alt);  %Continually at logical 1s for when we are at the altitudes where we exepct to be at waypoints
elseif size(waypt_alt,2)==2 %List of altitude bounds
    for i = 1:size(waypt_alt,1)
        msk_at_alt = msk_at_alt|...
            (alt>=waypt_alt(i,1)&...
             alt<=waypt_alt(i,2));  %Continually at logical 1s for when we are at the altitudes where we exepct to be at waypoints
    end
else
    error('UAVRT:waypointaltdef','Altitude definition was not a single scalar nor an nx2 matrix. ')
end

%Create initial waypoint candidate index list and mask. Note that these
%include points where the dwell times might be outside those specified by
%the user. 
    msk_at_waypt = [(V_smooth<v_thresh_at_wypt &...
                        time(1:end-1)>t_start & ...
                        time(1:end-1)<t_end&...
                        msk_at_alt(1:end-1));...
                        0];%Tack on a zero to make sure it is the same length as the original data record


    
    %Define entering waypoint when diff of mask is positive and leaving
    %when diff of mask is negative
    waypt_in_logic = diff(msk_at_waypt)>0; 
    waypt_out_logic = diff(msk_at_waypt)<0; 
    waypt_in_ind = find(waypt_in_logic);
    waypt_out_ind = find(waypt_out_logic);
    waypt_in_time = time(waypt_in_logic);
    waypt_out_time = time(waypt_out_logic);
    %Find point where there dwell time wasn't right.... (like take off and
    %landing)
	waypt_dt = waypt_out_time-waypt_in_time;
    
% Create the time limits for the autodetected waypoint dwell times. If
% waypt_dwell is a numeric vector, this will be passed. 
    if strcmp(waypt_dwell,'auto')     
    median_dwell = median(waypt_dt);
    std_dwell = std(waypt_dt);
    %Median +/- 3x the standard deviation;
    %Use the median to reduce likelyhood of grabbing an outlier
    waypt_dwell = median_dwell+3*[-std_dwell std_dwell];
    end

    waypt_clip_logic = waypt_dt>=min(waypt_dwell)&waypt_dt<=max(waypt_dwell);
    %And then remove them from the list of waypoints
    waypt_in_ind = waypt_in_ind(waypt_clip_logic);
    waypt_out_ind = waypt_out_ind(waypt_clip_logic); 
    waypt_in_time = waypt_in_time(waypt_clip_logic);
    waypt_out_time = waypt_out_time(waypt_clip_logic); 
    
    
% Create a plot of the velocity, altitude, and waypoint selects if
% requested
if strcmp(plot_control,'plot')     
     %figure; 
     %subplot(2,1,1)
     %axes(ax2plot)
     plot(ax2plot(1),time(1:end-1),V_smooth/max(V_smooth));hold(ax2plot(1),'on'); 
     plot(ax2plot(1),time,alt/max(alt));
     legend_strings{1} = ['Normalized smoothed vel. (max=',num2str(max(V_smooth)),'m/s)'];
     legend_strings{2} = ['Normalized smoothed alt. (max=',num2str(max(alt)),'m/s)'];
     waypt_hand = zeros(1,length(waypt_in_ind));
     for i = 1:length(waypt_in_ind)
         %save the handles so we can get the colors later to stay
         %consistent. 
         waypt_hand(i) = plot(ax2plot(1),time(waypt_in_ind(i):waypt_out_ind(i)),zeros(length(waypt_in_ind(i):waypt_out_ind(i)),1),'.');hold(ax2plot(1),'on');
         legend_strings{i+2} = ['Waypt #',num2str(i),' here'];
     end
     xlabel(ax2plot(1),'Time (s)')
     legend(ax2plot(1),legend_strings)
end

%This block create the waypoint summary list matrix  
msk_at_waypt_2 = zeros(size(time));%a new msk that exclude the points that were originally classified as waypoints but didn't have the right dwell time. 
waypt = zeros(length(waypt_in_ind),5); %preallocate array
    for i = 1:length(waypt_in_ind) % for each waypoint
       ind_range =  waypt_in_ind(i):waypt_out_ind(i);
       waypt(i,1) = mean(x_pos(ind_range));
       waypt(i,2) = mean(y_pos(ind_range));
       waypt(i,3) = mean(alt(ind_range));
       waypt(i,4) = waypt_in_time(i);
       waypt(i,5) = waypt_out_time(i);
       
       msk_at_waypt_2(ind_range) = 1;
    end
    
  
a_sum(:,10) = msk_at_waypt;  %send out the logic of being at a waypt
a_sum(:,11) = msk_at_alt;    %send out the logic of being at alt. This may be useful for transect flights. 

%Plot the 3d position of the flight and waypoints, along with a 2d ground
%projection.
if  strcmp(plot_control,'plot')      
    %figure
    %subplot(2,1,2) %plot the x-y trac
    plot3(ax2plot(2),a_sum(:,8),a_sum(:,9), a_sum(:,5),'Color',[0 0 1]);hold(ax2plot(2),'on')
    plot(ax2plot(2),a_sum(:,8),a_sum(:,9),'Color',[0.5 0.5 1]); 
    
    %Plot the 
    plot3(ax2plot(2),a_sum(find(msk_at_waypt_2),8),a_sum(find(msk_at_waypt_2),9),a_sum(find(msk_at_waypt_2),5),'o','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
    plot(ax2plot(2),a_sum(find(msk_at_waypt_2),8),a_sum(find(msk_at_waypt_2),9),'o','MarkerEdgeColor',[1 0.5 0.5],'MarkerFaceColor',[1 0.5 0.5]);
    
    the_waypt_colors = zeros(length(waypt_in_ind),3);%Preallocate array of the colors previous used for waypoint plotting
    the_waypt_colors_lighter = the_waypt_colors;
    for i = 1:length(waypt_in_ind)
        the_waypt_colors(i,:) = get(waypt_hand(i),'Color');
        the_waypt_colors_lighter(i,:) = the_waypt_colors(i,:) + ([1 1 1]-the_waypt_colors(i,:))*0.5; %Make 50% closer to white, in the same tone. 
        plot3(ax2plot(2),a_sum(waypt_in_ind(i):waypt_out_ind(i),8),a_sum(waypt_in_ind(i):waypt_out_ind(i),9),a_sum(waypt_in_ind(i):waypt_out_ind(i),5),'o','MarkerEdgeColor',the_waypt_colors(i,:),'MarkerFaceColor',the_waypt_colors(i,:));
        plot(ax2plot(2),a_sum(waypt_in_ind(i):waypt_out_ind(i),8),a_sum(waypt_in_ind(i):waypt_out_ind(i),9),'o','MarkerEdgeColor',the_waypt_colors_lighter(i,:),'MarkerFaceColor',the_waypt_colors_lighter(i,:));
     end

    xlabel(ax2plot(2),'X position North+ from home (m)')
    ylabel(ax2plot(2),'Y position East+ from home (m)')
    zlabel(ax2plot(2),'Alt. AGL from home (m)')
    %set(gca,'View',[90 90],'xdir','reverse')
    set(ax2plot(2),'View', [75 30],'xdir','reverse')
    axis(ax2plot(2),'equal')
    grid(ax2plot(2),'on')
end

 home = 180/pi*[lat_home, lon_home];%In degrees
end