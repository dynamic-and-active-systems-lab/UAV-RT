from datetime import datetime
import os

##########################################################################
# Variables
cwd = os.getcwd()  # Gets current directory, allows any user name to specified at setup
log_path_prefix = cwd[:-4] + '/CURRENT_DATA/'


##########################################################################
# Function

def gnu_output(proc, sock_out, UDP_IP_out, UDP_PORT_out):
    log_file_name = log_path_prefix + "log - " + datetime.now().strftime("%Y-%m-%d-T%H_%M_%S") + ".txt"
    log = open(log_file_name, "a", 0)
    while True:
        message = proc.stdout.readline()
        log.write(message)  # Saves to local file
        sock_out.sendto("//TRM:" + message, (UDP_IP_out, UDP_PORT_out))  # Sends Message to MATLAB
        print message,  # Prints message to terminal
