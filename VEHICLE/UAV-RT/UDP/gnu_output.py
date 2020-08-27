from datetime import datetime
import time
import os

##########################################################################
# Variables
cwd = os.getcwd()  # Gets current directory, allows any user name to specified at setup
log_path_prefix = cwd + '/UAV-RT/CURRENT_DATA/'

##########################################################################
# Function

def gnu_output(proc, sock_out, UDP_IP_out, UDP_PORT_out):
    log_file_name = log_path_prefix + "log - " + datetime.now().strftime("%Y-%m-%d-T%H_%M_%S") + ".txt"
    log = open(log_file_name, "a", 0)
    os.nice(10)  #Decrease the priority of this task
    while True:
        message = str(time.time())+': '+proc.stdout.readline()
        log.write(message)  # Saves to local file
        #There is a known issues when a ground control station is connected to the
        #Pixhawk as well as dronekit. Dronekit takes issue with hearbeats from the
        #ground control station and the pixhawk and throw the exceptions that we
        #silently ignore below. We have never seen these cause any problems. 
        if not ">>> Exception in message handler for HEARTBEAT" in message and \
        not ">>> mode 0 not available on mavlink definition" in message:
            sock_out.sendto("//TRM:" + message, (UDP_IP_out, UDP_PORT_out))  # Sends Message to MATLAB
            print message,  # Prints message to terminal
    
    '''
        #Below is where we check sucessfule connection and for known errors and pass them up to the parent process so that they can be properly dealt with.
        if "Using AirSpy NOS" in message:
            conn_to_parent.send('AIRSPY RUNNING')#send the message to the parent process
        if "Failed to stop RX streaming (-1000)" in message:
            conn_to_parent.send(message)#send the message to the parent process
        if "FATAL: Failed to open AirSpy device" in message:
            conn_to_parent.send('COULD NOT CONNECT TO AIRSPY - CHECK CONNECTION')#send the message to the parent process
    '''

