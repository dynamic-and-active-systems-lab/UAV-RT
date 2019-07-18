from Logger import Logger, Logger2
import sys
import socket
from multiprocessing import Process
import subprocess
import time
from datetime import datetime
from gnu_output import gnu_output
import os
import shutil
from config_read import config_read

time.sleep(30)  # delay allows network connections to be initialized on startup
#############################################################################
# VARIABLES
UDP_IP_in = "10.42.0.1"  # UDOO IP address
# UDP_IP_in = "134.114.64.117"  # UDOO IP address
UDP_PORT_in = 9090
UDP_PORT_out = 9091
cwd = os.getcwd()  # Gets current directory, allows any user name to specified when setting up UDOO
log_path_username = cwd  # directory should have a similar format to:/home/USERNAME
gnu_path = log_path_username + '/UAV-RT/GNU_RADIO/top_block.py'  # location of GNU Radio Script
path_prefix = log_path_username + '/UAV-RT/FLIGHT_DATA/'  # flight data final location path
current_data = log_path_username + '/UAV-RT/CURRENT_DATA/'  # file path for the latest flight data
config_path = log_path_username + '/UAV-RT/curr_config.uavrt'  # file path to configuration file

#############################################################################
# Intializes UDP Receiving IP Port
sock = socket.socket(socket.AF_INET,  # Internet
                     socket.SOCK_DGRAM)  # UDP
sock.bind((UDP_IP_in, UDP_PORT_in))

# Saves script print output for debugging/future reference
sys.stdout = Logger()
sys.stderr = Logger2()

##########################################################################
# MAIN PROGRAM
print 'Initializing UDP Connection Server'
while True:
    data, addr = sock.recvfrom(1024)  # buffer size is 1024 bytes

    ########################################
    # Allows user to send IP address to setup output port
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

    ########################################
    # Allows user to start GNU Radio
    if '//CMD:start' in data:
        # Info Command lets MATLAB know it received the start command
        sock_out.sendto('//INFO:start\n', (UDP_IP_out, UDP_PORT_out))

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

        # Start flight and process to send data over UDP
        sub_p = subprocess.Popen(['python', gnu_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        p = Process(target=gnu_output, args=(sub_p, sock_out, UDP_IP_out, UDP_PORT_out))
        p.start()

    ########################################
    # Allows user to stop GNU Radio
    if '//CMD:stop' in data:
        # Info Command lets MATLAB know it received the start command
        sock_out.sendto('//INFO:stop\n', (UDP_IP_out, UDP_PORT_out))

        print('Ending Program')  # print terminal output

        # terminate GNU Radio and process sending GNU Radio output to UDP
        sub_p.terminate()
        p.terminate()
        sock_out.sendto('//TRM:GNU Radio Ended\n', (UDP_IP_out, UDP_PORT_out))

        # Move files for flight
        files = os.listdir(current_data)
        for f in files:
            shutil.move(current_data+f, flight_folder_name)

        # Inform User the data is succesfully moved
        print 'Data Succesfully Moved'
        sock_out.sendto('//TRM:Data Succesfully Moved\n', (UDP_IP_out, UDP_PORT_out))

    ########################################
    # Allows user to determine if connection is alive
    if '//CMD:ECHO' in data:
        sock_out.sendto('//INFO:ECHO\n', (UDP_IP_out, UDP_PORT_out))
