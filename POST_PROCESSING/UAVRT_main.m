%% UAV-RT_main
% Version 1.0
% 2019-02-19
%
% UAV-RT_main processes flight and radio data from the UAV-RT system and
%displays results in a Matlab figure window, a .kml file for GoogleEarth, 
%and a text output file. The program assumes that flight logs and radio 
%data begin at the same time. Note thatt n the current system their 
%timestamps in the filenames may look different, but these are the times 
%when the files are initialized, not when the recordings starts. The 
%writing to those files does begin at the same time. 


%Author: Michael Shafer
%%************************************************************************

%%Input parameter:
%flt_data_flnm          Location of the flight data record
%sdr_data_flnm          Location of the IQ radio data record
%v_thresh_at_wypt       Velocity below which the vehicle should be at a
%                       waypoint
%Fs                     Sample rate of the IQ data
%pulse_dur              Expected pulse duration (s)
%pulse_rep              Expected pulse repetition rate (s)
%filt_band              Width of filter (Hz)
%savestring             a char array that will begin each output filename
%results_location       place where you want the results saved. 


%OUTPUTS
%savestring.txt         Summary files of processed results
%savestring.kml         Summary kml file

clear
clc

results_location = 'OUTPUT_FILES';

%Control variables
v_thresh_at_wypt = 0.5;
Fs = 48000;
pulse_dur = 0.02;
pulse_rep = 1.3;
filt_band = 160;


%% Example data for 61m altitude flight along 5 waypoints. 
flt_data_flnm = 'Example_data/flt_data_LINE-61m.txt';
sdr_data_flnm = 'Example_data/sdr_data_LINE-61m.dat';
waypt_dwell_time = 'auto';
savestring = 'EXAMP-61m';


%% FLIGHT DATA PROCESSING
%This was the already processed mat file: '/Users/mws22/Google Drive/DASL Drive/RESEARCH_AREAS/UAV/Testing Data - Documentation/2018-03-02-Lake Mormon range test with R-2 and circle test_using GNU sync/flight1_T10.17.27.800.mat';
% flt_data_flnm = '/Users/mws22/Google Drive/DASL Drive/RESEARCH_AREAS/UAV/Testing Data - Documentation/2018-03-02-Lake Mormon range test with R-2 and circle test_using GNU sync/Raw Data/test_data_telem-2018-03-02_10_17_27.txt';

%To consider only a portion of the record, you only need to change the time
%selection in this function. 
[veh_states,waypt,latlonhome] = flightprep(flt_data_flnm,'all','noplot','waypt_time',waypt_dwell_time);%'plot',1 'waypt_time',[80 87],[40 60]
waypt_time = waypt(:,end-1:end);%Grab the last two columns to create a more readable code with waypt_times [start end] for each waypoint (in the rows)
num_of_waypts = size(waypt,1);

%%  RADIO DATA INPUT AND PROCESSING

%RADIO DATA PREPARATION AND PULSE DETERMINATION
radio_out = radioprep(sdr_data_flnm,Fs,pulse_dur,'plot','none','fmax', 2000);%,'spectro',1);
    f = abs(radio_out{1}); %This is the frequency where pulses were detected.
    SDR_raw = radio_out{2};
%Filter and pull out maximums
data_out = cleandata(Fs,pulse_dur,f,filt_band,SDR_raw);
    t_g      = data_out{1};
    data_abs = abs(data_out{2});
%Find pulses within signal amplitude data
pulse_data = pulsefind(1/pulse_dur,pulse_dur,pulse_rep,data_abs);
    %Sometimes pulses are found when we don't have flight data (after
    %landing). Eliminate them with the mask that checks time stamps. 
    valid_pulse_msk = (pulse_data{1}>=veh_states(1,1)&pulse_data{1}<=veh_states(end,1));
    pulse_times = pulse_data{1}(valid_pulse_msk);
    pulse_amp   = pulse_data{2}(valid_pulse_msk);
    pulse_pow   = pulse_data{2}(valid_pulse_msk).^2;%power
	num_of_pulses = length(pulse_times);
%figure;plot(pulse_times,pulse_amp,'.','Markersize',10); hold on;plot(pulse_times,pulse_amp);xlabel('Time(s)');ylabel('Pulse amplitude');


%% FLIGHT DATA AND RADIO DATA FUSION

veh_out = vehstates(pulse_times, veh_states,waypt_time);%
    veh_out{1}(:,4) = wrapTo360(veh_out{1}(:,4)); %Wrap heading to 360. No longer need total angle.
    pulse_veh_states = veh_out{1};%List of vehilce states at the moment of pulse reception
    pulse_x   = veh_out{1}(:,8); % X position when each pulse is received
    pulse_y   = veh_out{1}(:,9); % Y position when each pulse is received
    pulse_yaw = veh_out{1}(:,4); % Yaw (heading in deg)  when each pulse is received
    pulse_alt = veh_out{1}(:,5); % Altitude when each pulse is received
    pulse_waypt_num = veh_out{2};% The waypoint the vehicle was at when this pulse is received. NaN if not at a waypoint

%This code is use to find the maximum and average pulse power and amp at
%each waypoint.
    waypt_max_pulse_pow = zeros(num_of_waypts,1);
    waypt_avg_pulse_pow = waypt_max_pulse_pow;
    waypt_max_pulse_amp = waypt_max_pulse_pow;
    waypt_avg_pulse_amp = waypt_max_pulse_pow;

for i = 1:num_of_waypts %This finds the maximum pulse signal amplitude at each waypoint
    waypt_max_pulse_pow(i) =  max(pulse_pow(pulse_waypt_num == i));
    waypt_avg_pulse_pow(i) = mean(pulse_pow(pulse_waypt_num == i));
    waypt_max_pulse_amp(i) =  max(pulse_amp(pulse_waypt_num == i));
    waypt_avg_pulse_amp(i) = mean(pulse_amp(pulse_waypt_num == i));
end

if num_of_waypts~=0
    %% DOA ESTIMATION
    %                (Pulse Amplitudes, Pulse Yaws, Pulse waypoint number, power or amp?, scaling,plot control)
    [doaout] = doapca(pulse_amp,pulse_yaw,pulse_waypt_num,'power','linear','noplot');%All inputs required
        DOA_tau = doaout{2};                                             %Tau values of each bearing
        DOA_calc_deg_N_CW = doaout{1};                                   %DOA degree angle from N with postiive CW -i.e. a regular compas bearing.
        DOA_calc_rad_N_CW = DOA_calc_deg_N_CW*pi/180;                    %DOA radian angle from N with postiive CW -i.e. a regular compas bearing.
        DOA_calc_rad_E_CCW = wrapToPi(pi/2-wrapToPi(DOA_calc_rad_N_CW)); %DOA radian angle from E with postiive CCW -i.e. standard math angle definition
        DOA_calc_deg_E_CCW = 180*pi*DOA_calc_rad_E_CCW;                  %DOA degree angle from E with postiive CCW -i.e. standard math angle definition
    
    %% LOCALIZATION
    %The localization algorithms use East as X and N as Y, so we enter them in
    %that way here and name the outputs yx to remind use that the entries are
    %flipped. We'll switch back to [x y] after some conversions.
    [pos_yNxE_CM,intermap_yNxE] = localize_cm(waypt(:,2),waypt(:,1),DOA_calc_deg_N_CW,waypt_max_pulse_pow, DOA_tau,0.35,30);
    [pos_yNxE_MLE]              = localize(waypt(:,2),waypt(:,1),DOA_calc_deg_N_CW,'MLE');
    [pos_yNxE_RMR]              = localize(waypt(:,2),waypt(:,1),DOA_calc_deg_N_CW,'RMR');
    [pos_yNxE_MEST]             = localize(waypt(:,2),waypt(:,1),DOA_calc_deg_N_CW,'MEST');
    [pos_yNxE_MEDIAN]           = median([pos_yNxE_CM,pos_yNxE_MLE,pos_yNxE_RMR,pos_yNxE_MEST],2);%Median of all. Not mean in case there is an outlier.
     
    %% CONVERT POSITIONS
        %Localization results in lat/lons
        latlon_CM     = xy2latlon(latlonhome,pos_yNxE_CM');
        latlon_MLE    = xy2latlon(latlonhome,pos_yNxE_MLE');
        latlon_RMR    = xy2latlon(latlonhome,pos_yNxE_RMR');
        latlon_MEST   = xy2latlon(latlonhome,pos_yNxE_MEST');
        latlon_MEDIAN   = xy2latlon(latlonhome,pos_yNxE_MEDIAN');
        %Localization results in cylindrical coordinates (radial distance
        %and bearing)
        cyl_CM     = [norm(pos_yNxE_CM),wrapTo360(180/pi*atan2(pos_yNxE_CM(1),pos_yNxE_CM(2)))];
        cyl_MLE    = [norm(pos_yNxE_MLE),wrapTo360(180/pi*atan2(pos_yNxE_MLE(1),pos_yNxE_MLE(2)))];
        cyl_RMR    = [norm(pos_yNxE_RMR),wrapTo360(180/pi*atan2(pos_yNxE_RMR(1),pos_yNxE_RMR(2)))];
        cyl_MEST   = [norm(pos_yNxE_MEST),wrapTo360(180/pi*atan2(pos_yNxE_MEST(1),pos_yNxE_MEST(2)))];
        cyl_MEDIAN   = [norm(pos_yNxE_MEDIAN),wrapTo360(180/pi*atan2(pos_yNxE_MEDIAN(1),pos_yNxE_MEDIAN(2)))];
        %Waypoints in lat/lon
        %Switch the x-y because the lat/lon converter uses xEyN
        waypt_latlon = xy2latlon(latlonhome,waypt(:,2:-1:1));
        pulse_latlon = xy2latlon(latlonhome,[pulse_y,pulse_x]);
        pulse_cyl   = [hypot(pulse_y,pulse_x),wrapTo360(180/pi*atan2(pulse_y,pulse_x))];
        %Flip the x and y's to get back to X-North Y-East frame.
        pos_xNyE_CM   = flip(pos_yNxE_CM);intermap_xNyE = intermap_yNxE(:,[2,1,3:end]);
        pos_xNyE_MLE  = flip(pos_yNxE_MLE);
        pos_xNyE_RMR  = flip(pos_yNxE_RMR);
        pos_xNyE_MEST = flip(pos_yNxE_MEST);
        pos_xNyE_MEDIAN = flip(pos_yNxE_MEDIAN);
    
else
    warning('No waypoints detected. Skipping DOA and Localization')
end

%% MAKE THE RESULTS MAP

figure

%Plot entire flight as a blue line
plot3(veh_states(:,8),veh_states(:,9), veh_states(:,5),'Color',[0 0 1]); hold on %3D points
plot(veh_states(:,8),veh_states(:,9),'Color',[0.9 0.9 1]); hold on            %on the ground

%Plot pulses as green circles with radii proportional to received
%amplitude
for i = 1:length(pulse_x)
    plot3(pulse_x(i),pulse_y(i),pulse_alt(i),'o','Markersize',30*pulse_amp(i)/max(pulse_amp),'MarkerEdgeColor',[0.5 1 0.5]); hold on;
end

if num_of_waypts~=0
    flight_delta_dist = sqrt(range(veh_states(:,8)).^2+range(veh_states(:,9)).^2);
    vect_lengths =10*flight_delta_dist*waypt_max_pulse_pow/max(waypt_max_pulse_pow);
    
    %Plot waypts flight as a red dots
    plot3(waypt(:,1),waypt(:,2),waypt(:,3),'o','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
    plot(waypt(:,1),waypt(:,2),'o','MarkerEdgeColor',[1 0.9 0.9],'MarkerFaceColor',[1 0.9 0.9]);
    quiver3(waypt(:,1),waypt(:,2),waypt(:,3),vect_lengths.*cos(DOA_calc_rad_N_CW),vect_lengths.*sin(DOA_calc_rad_N_CW),zeros(size(DOA_calc_rad_N_CW)),'Color',[0.5 0.5 1])
    quiver(waypt(:,1),waypt(:,2),vect_lengths.*cos(DOA_calc_rad_N_CW),vect_lengths.*sin(DOA_calc_rad_N_CW),'Color',[0.9 0.9 1])
    %Plot the localization results
    plot(pos_xNyE_CM(1),pos_xNyE_CM(2),'bs','Markersize',10)
    plot(pos_xNyE_MLE(1),pos_xNyE_MLE(2),'b+','Markersize',10)
    plot(pos_xNyE_RMR(1),pos_xNyE_RMR(2),'bo','Markersize',10)
    plot(pos_xNyE_MEST(1),pos_xNyE_MEST(2),'b*','Markersize',10)
end

xlabel('X position +North from home (m)')
ylabel('Y position +East from home (m)')
zlabel('Alt. AGL from home (m)')
axis equal
grid on

set(gca,'View',[90 90],'xdir','reverse')


%% KML Generation

temp_kml_dir_name = 'temp_KML';
 mkdir(temp_kml_dir_name)
 
 %Write the flight path to a kml file and extract the needed string data
 kmlwriteline('temp_KML/flight_path.kml',veh_states(:,6),veh_states(:,7),veh_states(:,5),'AltitudeMode','relativeToGround','Color',[0.5 0.5 1],'Name','Flight Path','Alpha',1);
 flight_path_text = kmlextract('temp_KML/flight_path.kml'); 
 %Write the pulses to a kml file and encapsulate in a folder
 for i = 1:num_of_pulses;pulse_blank_labels{i} = ' ';end
 pulse_icon_size = 3*pulse_pow/max(pulse_pow); %Range is 2*(0-1);
 parula_10 = jet(10);
 pulse_colors = parula_10(round(1+9*round(pulse_pow/max(pulse_pow),1)),:);
 kmlwritepoint('temp_KML/pulse_pow.kml',pulse_latlon(:,1),pulse_latlon(:,2),pulse_alt,'Name',pulse_blank_labels,'Color',pulse_colors,'Icon','kml_files/wht-glo.png','Iconscale',pulse_icon_size,'AltitudeMode','relativeToGround');
 %Encapsulate pulses into a folder for the eventual kml file
 folder_list_pulse = {'temp_KML/pulse_pow.kml'};
 folder_text_pulse = kmlfolder(folder_list_pulse,'Received Pulses');

 %Only execute this if there were waypoints the vehicle flew to
 if num_of_waypts~=0    
     %Waypoints
     for i = 1:num_of_waypts;waypt_name_list{i} = sprintf('WP%i',i);end
     kmlwritepoint('temp_KML/vehicle_waypts.kml',waypt_latlon(:,1),waypt_latlon(:,2),0*waypt(:,3),'AltitudeMode','relativeToGround','Icon','kml_files/logo_blue.png','Name',waypt_name_list,'IconScale',1*ones(1,num_of_waypts));
     %Location Estimates
     kmlwritepoint('temp_KML/est_pos_cm.kml',latlon_CM(1),latlon_CM(2),'Name','CM Est.','Color',[0.5 0.5 1],'Icon','kml_files/wht-CM.png');
     kmlwritepoint('temp_KML/est_pos_mle.kml',latlon_MLE(1),latlon_MLE(2),'Name','MLE Est.','Color',[0.5 0.5 1],'Icon','kml_files/wht-MLE.png');
     kmlwritepoint('temp_KML/est_pos_rmr.kml',latlon_RMR(1),latlon_RMR(2),'Name','RMR Est.','Color',[0.5 0.5 1],'Icon','kml_files/wht-RMR.png');
     kmlwritepoint('temp_KML/est_pos_mest.kml',latlon_MEST(1),latlon_MEST(2),'Name','M-Est.','Color',[0.5 0.5 1],'Icon','kml_files/wht-MEST.png');
     kmlwritepoint('temp_KML/est_pos_med.kml',latlon_MEDIAN(1),latlon_MEDIAN(2),'Name','Median Est.','Color',[0.5 0.5 1],'Icon','kml_files/wht-MED.png');
     %Bearing Estimates
     %Each draw distance is the distance from that waypoint to the estimated tag location 
     draw_dists = sqrt((pos_xNyE_MEDIAN(1)-waypt(:,1)).^2+(pos_xNyE_MEDIAN(2)-waypt(:,2)).^2); 
     %Plot at alt of 1 m.
     kmlbearing('temp_KML/bear_lines.kml',waypt_latlon(:,1),waypt_latlon(:,2),ones(size(waypt(:,3))),DOA_calc_deg_N_CW,draw_dists,'blue',1);
     
     
     %Encapsulate each group into folders for the eventual kml file
     folder_list_est_pos = {'temp_KML/est_pos_cm.kml','temp_KML/est_pos_mle.kml','temp_KML/est_pos_rmr.kml','temp_KML/est_pos_mest.kml','temp_KML/est_pos_med.kml'};
     folder_text_est = kmlfolder(folder_list_est_pos,'Location Estimates'); 
     folder_list_waypts = {'temp_KML/vehicle_waypts.kml'};
     folder_text_waypts = kmlfolder(folder_list_waypts,'Waypoints');
     folder_list_bear = {'temp_KML/bear_lines.kml'};
     folder_text_bear = kmlfolder(folder_list_bear,'Bearings');
 end

doc_text_pre = {flight_path_text,folder_text_est,folder_text_waypts,folder_text_bear,folder_text_pulse};
doc_text_compile = kmlwriteprep(doc_text_pre,savestring);

fid = fopen([results_location,'/',savestring,'-MAP.kml'],'wt');
fprintf(fid, doc_text_compile);
fclose(fid);
rmdir(temp_kml_dir_name,'s') %Remove the temporary folder and it contents



%% PREPARE SUMMARY FILE FROM RESULTS
if num_of_waypts~=0
    %Waypoint matrix will include: waypt #, DOA from N, waypoint relative pulse power, waypoint lat, waypoint lon, waypoint x, waypoint y, waypoint altitude,  time start waypoint, time end waypoint
    waypt_out_mat = [(1:num_of_waypts)',DOA_calc_deg_N_CW,waypt_max_pulse_pow/max(waypt_max_pulse_pow),waypt_latlon,waypt];
    
    %Localiazation matrix will in clude each methods: est. dist to tag,  bearing to tag, est-lat est-lon, est-x, est-y, est-lat est-lon
    loc_out_mat_names = ['MEDIAN';'    CM';'   MLE';'   RMR';'  MEST'];
    loc_out_mat_data = [cyl_MEDIAN,  latlon_MEDIAN,  pos_xNyE_MEDIAN';...
        cyl_CM,  latlon_CM,  pos_xNyE_CM';...
        cyl_MLE, latlon_MLE, pos_xNyE_MLE';...
        cyl_RMR,  latlon_RMR,  pos_xNyE_RMR';...
        cyl_MEST, latlon_MEST, pos_xNyE_MEST'];
    loc_out_mat = [cellstr(loc_out_mat_names), num2cell(loc_out_mat_data)].';
    
    waypt_header_text = ['#  DOA(deg-N)','\t','WP RSS','\t','Lat','\t \t','Lon','\t \t ','x(m)','\t \t','y(m)','\t \t','Alt(m)','\t','WP t_0','\t \t','WP t_end \r\n ------------------------------------------------------------------------------------------------------------------------\r\n'] ;
    loc_header_text = ['Method  Dist.(m) Deg-N','\t','Lat','\t \t','Lon','\t \t','x(m)','\t','y(m)','\r\n -----------------------------------------------------------------------\r\n'] ;
end

strongest_pulse_index = find(pulse_amp==max(pulse_amp(2:end-1)));

output_header =  ['Results of input files:','\n',...
    'Flight data: ',flt_data_flnm,'\n',...
    'Radio data: ',sdr_data_flnm,'\n',...
    '\n',...
    '\n',...
    '\n',...
    num2str(num_of_pulses),' radio pulses found at ',num2str(f),' Hz from center frequency','\n',...
    'Expected approximately ', num2str(round((veh_states(end,1) - veh_states(1,1))/pulse_rep)),' pulses for a ',num2str((veh_states(end,1) - veh_states(1,1))/60),' minute flight.','\n',...
    'Waypoints found: ',num2str(num_of_waypts),'\n',...
    'Strongest pulse found at [x,y] = [',num2str([pulse_x(strongest_pulse_index),pulse_y(strongest_pulse_index)]),']\n',...%ignore the first pulse becasue sometimes the first pulse is erroneous as the SDR is getting started up. The first second or so of data should typically be ignored.
    '\t \t \t [lat lon] = [',num2str([pulse_latlon(strongest_pulse_index,1),pulse_latlon(strongest_pulse_index,2)]),'] \n',...
    '\t \t \t [range bearing] = [',num2str([pulse_cyl(strongest_pulse_index,1),pulse_cyl(strongest_pulse_index,2)]),'] [meters, degrees-N]) \n'];

pulse_header_text = ['#','\t','Time(s)','\t','Pow-mW','\t \t','Lat','\t \t','Lon','\t \t ','x(m)','\t \t','y(m)','\t \t','Alt(m)',' \r\n ------------------------------------------------------------------------------------------------------------------------\r\n'] ;
pulse_out_mat = [(1:length(pulse_times))',pulse_times',pulse_pow',pulse_latlon(:,1),pulse_latlon(:,2),pulse_x,pulse_y,pulse_alt];


%Write the output file
fileID = fopen([results_location,'/',savestring,'-SUMMARY.txt'],'w');
fprintf(fileID,output_header);
if num_of_waypts~=0
fprintf(fileID,'LOCALIZATION RESULTS: \r\n \r\n');
fprintf(fileID,loc_header_text);
fprintf(fileID,'%s %6.1f \t %.1f \t %.8f \t %.8f \t %.1f \t %.1f \r\n',loc_out_mat{:});
fprintf(fileID,'\r\n \r\n \r\n');
fprintf(fileID,'WAYPOINT DATA LIST: \r\n \r\n');
fprintf(fileID,waypt_header_text);
fprintf(fileID,'%2u %6.2f \t %2.1f \t %.8f \t %.8f \t %7.1f \t %7.1f \t %.1f \t %7.2f \t %7.2f\r\n',waypt_out_mat');
end
fprintf(fileID,'\r\n \r\n \r\n');
fprintf(fileID,'RECEIVED PULSE DATA LIST: \r\n \r\n');
fprintf(fileID,pulse_header_text);
fprintf(fileID,'%5u %7.2f \t %4.2e \t %.8f \t %.8f \t %7.1f \t %7.1f \t %.1f\r\n',pulse_out_mat');
fclose(fileID);


    
    


