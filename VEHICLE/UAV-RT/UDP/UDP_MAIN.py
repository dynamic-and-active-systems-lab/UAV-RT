from Logger import Logger, Logger2
import sys
import socket
from multiprocessing import Process, Pipe
import subprocess
import time
from datetime import datetime
from gnu_output import gnu_output
import os
import shutil
from config_read import config_read
from threading import Thread

#########################
# This function will be used to start up the gnuradio script as a subprocess
def start_radio():
    # Create folder for data (note not moved into folder until after flight)
    flight_folder_name = path_prefix + "UAV-RT - " + datetime.now().strftime("%Y-%m-%d-T%H_%M_%S")
    os.makedirs(flight_folder_name)

    # Read Config file
    config_read(config_path, gnu_path, flight_folder_name)
    print('Reading Configuration File')
    sock_out.sendto('//TRM:Reading Configuration File\n', (UDP_IP_out, UDP_PORT_out))

    # Print output to user terminal and UDP output
    print('Starting GNU Radio - ') + time.strftime('%s')  # time stamps for reference if debugging
    sock_out.sendto('//TRM:Starting GNU Radio\n', (UDP_IP_out, UDP_PORT_out))
    
    #parent_conn, child_conn = Pipe()   #create pipe connection points for process
    try:        
        # Start flight and process to send data over UDP
        sub_p = subprocess.Popen(['python', gnu_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print 'Started subprocess'
        p = Process(target=gnu_output, args=(sub_p, sock_out, UDP_IP_out, UDP_PORT_out))
        print 'Created gnu_output'
        p.start()
        print 'Started gnu_output process'

    except:
        print('GNU Radio Script Error') + time.strftime('%s')  # time stamps for reference if debugging
        sock_out.sendto('//TRM:RADIO START ERROR!!!\n', (UDP_IP_out, UDP_PORT_out))
        stuff_to_return = [False,False,False,False]
    else:
        # Print output to user terminal and UDP output
         print('GNU Radio Script Running') + time.strftime('%s')  # time stamps for reference if debugging
         sock_out.sendto('//TRM:Radio Process Initiated\n', (UDP_IP_out, UDP_PORT_out))
         stuff_to_return = [p,sub_p,flight_folder_name]
    return stuff_to_return
    
time.sleep(10)  # delay allows network connections to be initialized on startup
                # may want to reduce this or comment out if connected to a 
                # wired network for testing
#############################################################################
# VARIABLES
UDP_IP_in = "10.42.0.1"  # UDOO IP address
# UDP_IP_in = "134.114.64.117"  # UDOO IP address
UDP_PORT_in = 9090
UDP_PORT_out = 9091
cwd = os.getcwd()  # Gets current directory, allows any user name to specified when setting up UDOO
log_path_username = cwd  # directory should have a similar format to:/home/USERNAME
gnu_path = log_path_username + '/UAV-RT/VEHICLE/UAV-RT/GNU_RADIO/top_block.py'  # location of GNU Radio Script
path_prefix = log_path_username + '/UAV-RT/VEHICLE/UAV-RT/FLIGHT_DATA/'  # flight data final location path
current_data = log_path_username + '/UAV-RT/VEHICLE/UAV-RT/CURRENT_DATA/'  # file path for the latest flight data
config_path = log_path_username + '/UAV-RT/VEHICLE/UAV-RT/curr_config.uavrt'  # file path to configuration file

#############################################################################
# Intializes UDP Receiving IP Port
sock = socket.socket(socket.AF_INET,  # Internet
                     socket.SOCK_DGRAM)  # UDP
sock.bind((UDP_IP_in, UDP_PORT_in))
#sock.setblocking(0) #Set to non-blocking so that loop can proceed if no data is available
#sock.settimeout(1)#Set a 1 second timeout on the receipt of data

# Saves script print output for debugging/future reference
sys.stdout = Logger()
sys.stderr = Logger2()

##########################################################################
# Create hearbeat background task
def heartbeat_task():
    os.nice(10)         #Decrease the priority of this task
    while not heartbeat_task_cancel:
        print 'HB '
        if 'parent_conn' in globals():
            gnu_message = parent_conn.recv()
            if "AIRSPY RUNNING" in gnu_message:
                sock_out.sendto('//STATUS: AIRSPY IS RUNNING'+'\n', (UDP_IP_out, UDP_PORT_out))
        if 'sock_out' in globals():  #Send if connection is active
            #print 'HB'
            #sock_out.sendto('//HB:'+str(time.time())+'\n', (UDP_IP_out, UDP_PORT_out))
            sock_out.sendto('//HB:\n', (UDP_IP_out, UDP_PORT_out))
            sock_out.sendto('//STATUS:'+status+'\n', (UDP_IP_out, UDP_PORT_out))
        time.sleep(1)
        
##########################################################################
# MAIN PROGRAM
print 'Initializing UDP Connection Server'
status = 'not running'
while True:
    #try
    data, addr = sock.recvfrom(1024)  # buffer size is 1024 bytes
    '''
    except:  #Send status on timeout
        print 'HB'            
        if 'sock_out' in locals():  #Send if connection is active
            sock_out.sendto('//HB:'+str(time.time())+'\n', (UDP_IP_out, UDP_PORT_out))
            sock_out.sendto('//STATUS:'+status+'\n', (UDP_IP_out, UDP_PORT_out))

        #t_next_HB = int(time.time())+HB_interval #reset the current integer time
        #print t_next_HB
    '''
    #else:
    ########################################
    # Allows user to send IP address to setup output port
    #Don't need to have sock_out defined for this, as this code will set up sock_out
    if '//IP' in data:
        start_point = data.find(':') + 1
        IP = data[start_point:-2]
        print 'Initializing UDP Output on', IP
        # UDP Sending IP and Port
        UDP_IP_out = IP
        sock_out = socket.socket(socket.AF_INET,  # Internet
                                 socket.SOCK_DGRAM)  # UDP
        # Info Command lets MATLAB know it received the IP command
        sock_out.sendto('//INFO:IP\n', (UDP_IP_out, UDP_PORT_out))
        # Tell MATLAB that connection was succesfully established
        sock_out.sendto('//TRM:Connection Succesful\n', (UDP_IP_out, UDP_PORT_out))
        #Setup and initiate the heart beat
        heartbeat_task_cancel=False
        HB_thread = Thread(target=heartbeat_task)
        HB_thread.start()

    ########################################
    # Allows user disconnect and properly shut down connections on this side.
    if '//CMD:disconnect' in data:
        print 'Disconnecting UDP Output'
        if 'sock_out' in locals():         
            # sock_out.shutdown(socket.SHUT_RDWR)
            sock_out.close()
            del sock_out   #Delete variable so we don't try to write to it again
            heartbeat_task_cancel=True #Kill the heartbeats
    ########################################
    # Must be connected with a valid output port to start or stop
    if 'sock_out' in locals():  #Only check for commands if we have first setup the output connection
        ########################################
        # Allows user to start GNU Radio
        if '//CMD:start' in data:
            # Info Command lets MATLAB know it received the start command
            sock_out.sendto('//INFO:start\n', (UDP_IP_out, UDP_PORT_out))
          
            if 'p' in locals(): # first check to see if the subprocess exists 
                    if p.is_alive(): #now check if the process is alive
                        print('Radio already running')
                        sock_out.sendto('//TRM:Radio already running\n', (UDP_IP_out, UDP_PORT_out))
                    else:
                        try:
                            print 'Trying to start radio up'
                            start_out = start_radio()
                            print 'DONE!'
                        except:
                            print 'There was an error starting up the radio'
                            status = 'not running'                        
                        else:
                            print 'Unpacking the processes to create variable'
                            p         = start_out[0]
                            sub_p     = start_out[1]
                            flight_folder_name = start_out[2]
                            #parent_conn = start_out[3]
                            status = 'running'

            else:  #if it wasn't already created, go ahead and start
                try:
                    print 'Trying to start radio up - p not yet in locals'
                    start_out = start_radio()
                    print 'DONE - Trying to start radio up - p not yet in locals'
                except:
                    print 'There was an error starting up the radio'
                    status = 'not running'                        
                else:
                    print 'Unpacking the processes to create variable'
                    p         = start_out[0]
                    sub_p     = start_out[1]
                    flight_folder_name = start_out[2]
                    #parent_conn = start_out[3]
                    status = 'running'
                    print 'DONE- Unpacking the processes to create variable'
 
            
        ########################################
        # Allows user to stop GNU Radio
        if '//CMD:stop' in data:
            # Info Command lets MATLAB know it received the start command
            sock_out.sendto('//INFO:stop\n', (UDP_IP_out, UDP_PORT_out))
            if 'p' in locals(): # first check to see if the subprocess exists 
                if p.is_alive(): #now check if the process is alive
                    print('Ending Program')  # print terminal output
                            
                    # terminate GNU Radio and process sending GNU Radio output to UDP
                    sub_p.terminate()
                    p.terminate()
                    #del parent_conn #Delete the connection so that HB task doesn't continue to try to read
                    status = 'not running'                        
                    sock_out.sendto('//TRM:GNU Radio Ended\n', (UDP_IP_out, UDP_PORT_out))
                    
                    # Move files for flight
                    files = os.listdir(current_data)
                    for f in files:
                        shutil.move(current_data+f, flight_folder_name)
                        
                    # Inform User the data is succesfully moved
                    print 'Data Succesfully Moved'
                    sock_out.sendto('//TRM:Data Succesfully Moved\n', (UDP_IP_out, UDP_PORT_out))
            else:
                sock_out.sendto('//TRM:Radio was not running.\n', (UDP_IP_out, UDP_PORT_out))
         ########################################
        # Allows user to determine if connection is alive
        if '//CMD:ECHO' in data:
            sock_out.sendto('//INFO:ECHO\n', (UDP_IP_out, UDP_PORT_out))

 
