import shutil
from tempfile import mkstemp
import os


def config_read(config_path, gnu_path, flight_folder_name):
    #############################################################################
    # Open Config File
    with open(config_path, 'r') as f_conf:
        lines = f_conf.readlines()

    # Read Config file variables
    for line in lines:
        ########################################
        # CONDITIONS FOR THE CONFIG VARIABLES
        # RF Gain
        if 'RFG' in line:
            start = line.find('\t')
            RFG = line[start+1:]
        # IF Gain
        if 'IFG' in line:
            start = line.find('\t')
            IFG = line[start+1:]
        # Frequecny
        if 'FT' in line:
            start = line.find('\t')
            FT = line[start+1:]
        # Sampling Rate
        if 'FRS' in line:
            start = line.find('\t')
            FRS = line[start+1:]
        ########################################

    f_conf.close()  # close configuration file

    # Copy config file for flight
    new_config_path = flight_folder_name + '/config.uavrt'
    shutil.copyfile(config_path, new_config_path)

    #############################################################################
    # Write GNU Radio Variables
    fh, temp_path = mkstemp()  # makes temp path for new file

    with os.fdopen(fh, 'w') as f_new:
        with open(gnu_path) as f:
            lines = f.readlines()
            for line in lines:
                ########################################
                # CONDITIONS FOR GNU VARIABLES
                # RF Gain
                if 'gainRF = gainRF =' in line:
                    line = line.replace(line[31:], RFG)
                # IF Gain
                if 'gainIF = gainIF =' in line:
                    line = line.replace(line[31:], IFG)
                # Frequecny
                if 'freq = freq =' in line:
                    line = line.replace(line[27:], FT[:-2] + 'e6\r\n')  # add e6 for Mhz
                # Sampling Rate
                if 'audio_rate = audio_rate =' in line:
                    line = line.replace(line[39:], FRS)
                ########################################
                # Write the new line
                f_new.write(line)

    shutil.move(gnu_path, gnu_path + '.bak')
    shutil.move(temp_path, gnu_path)
