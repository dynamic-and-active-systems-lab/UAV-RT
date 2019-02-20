#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018 <+YOU OR YOUR COMPANY+>.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

import numpy
from gnuradio import gr

import time
import os
import datetime
from dronekit import connect, VehicleMode
import subprocess
import sys
import platform
import math
import serial


class pixstream_source(gr.sync_block):
    """
    docstring for block pixstream_source
    """

    def __init__(self, telem_file_path):
        gr.sync_block.__init__(self,
                               name="pixstream_source",
                               in_sig=None,
                               out_sig=[(numpy.float32, 7)])

        # Check the OS to see what connection string to use to talk to Pixhawk

        osName = platform.system()
        if(osName == 'Windows'):
            connection_string = "COM3"  # for connection via USB
        elif(osName == 'Darwin'):  # Darwin is the result for a Mac
            connection_string = "/dev/cu.usbmodem1"
            # The Dronekit site doesn't have a / in the front of the conenction
            # but it is necessary
        else:
            connection_string = "/dev/ttyACM0"  # for connection via USB

        print(osName)

        # connection_string = "127.0.0.1:14550"

        # Now connect to Pixhawk
        print("Connecting to Pixhawk")
        self.vehicle = connect(connection_string, wait_ready=True, rate=10)
        print("Connected!")

        # Create the telemetry log file and add the header
        self.now = datetime.datetime.now()
        self.telem_file_name = telem_file_path+"/test_data_telem-" + \
            self.now.strftime("%Y-%m-%d-T%H.%M.%S") + ' - ' + time.strftime('%s') + ".txt"
        file = open(self.telem_file_name, "w")
        file.write("Time_s"+"\t""Heading_deg"+"\t"+"Pitch_rad"+"\t"+"Roll_rad" +
                   "\t"+"Yaw_rad"+"\t"+"Lat"+"\t"+"Long"+"\t"+"Alt_m_AGL"+"\n")
        file.close

        # I was playing around with variable save rates below, but it doesn't
        # currently work----
        # Create the class variable save_rate, so we can dictate how often we
        # save data in the log_file
        # self.save_rate = save_rate

    def write_telem(self, the_data):
        # This function is used to write the data to the tab delimited text file
        file = open(self.telem_file_name, "a")
        # print str(the_data)
        # Time, Heading, Pitch, Roll, Yaw, Lat, Lon, Alt-AGL
        file.write(str(the_data[0, 0])+"\t"+str(the_data[0, 1])+"\t"+str(the_data[0, 2])+"\t"+str(the_data[0, 3]) +
                   "\t"+str(the_data[0, 4])+"\t"+str(the_data[0, 5])+"\t"+str(the_data[0, 6])+"\t"+str(the_data[0, 7])+"\n")
        file.close

    def work(self, input_items, output_items):
        '''This method is called by the gnuradio code much more frequently than
        the 4-10 Hz that we want to save data. To limit the data rate, the
        method checks to see if the current time has a zero in the tenths place.
        If it does, it will call the 'write_telem' method to write the current
        telemetry data to the log file. If it isn't yet time to write the data,
        the method will sleep for a hundreth of a second and check again. The
        functionality of this method relies on the fact that the rest of the
        gnuradio code is running much faster and calling this function more
        frequently than every 0.01 seconds.
        '''
        # Determine the the time when the work method is called
        starttime = time.time()
        # Find the next 10th of a second. That is when we want to write data
        nexttime = (math.floor(starttime*10)+1)/10
        # nexttime = (math.floor(starttime*self.save_rate)+1)/self.save_rate
        # print "Start Time is: ", starttime
        # print "Next  Time is: ", nexttime
        nowtime = starttime  # nowtime is now. Will be updated later
        while nowtime <= nexttime:
            # <+signal processing here+>
            time.sleep(0.01)
            # time.sleep(1/self.save_rate*1/10)
            nowtime = time.time()
        # Once we get to this point (out of the while loop), we are at the
        # of the second when we want to write data....

        # Set up and fill the data vector
        data = numpy.zeros((1, 8))
        data[0, 0] = time.time()
        data[0, 1] = self.vehicle.heading
        data[0, 2] = self.vehicle.attitude.pitch
        data[0, 3] = self.vehicle.attitude.roll
        data[0, 4] = self.vehicle.attitude.yaw
        data[0, 5] = self.vehicle.location.global_relative_frame.lat
        data[0, 6] = self.vehicle.location.global_relative_frame.lon
        data[0, 7] = self.vehicle.location.global_relative_frame.alt
        # Write the data
        self.write_telem(data)

        # Below is where we stream the data out to gnu radio
        '''Setup the out variable that gnuradio is expecting to be the size that
        is requested. In my tests to this point, the output_items length
        is on the order of 1500 elements. I haven't figured out how to limit
        which would have made it able to just send the 'data' variable. It looks
        like you can change the output buffer of a block written in C++, but now
        one written in python. We have to use python because that is the API we
        have for Dronekit.
        '''
        out = output_items[0]
        out[:, 0] = self.vehicle.heading
        out[:, 1] = self.vehicle.attitude.pitch
        out[:, 2] = self.vehicle.attitude.roll
        out[:, 3] = self.vehicle.attitude.yaw
        out[:, 4] = self.vehicle.location.global_relative_frame.lat
        out[:, 5] = self.vehicle.location.global_relative_frame.lon
        out[:, 6] = self.vehicle.location.global_relative_frame.alt

        return len(output_items[0])
