#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Top Block
# Generated: Tue Jul 16 14:12:21 2019
##################################################


from datetime import datetime
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import osmosdr
import time
import timesync
import os  # MODIFIED


class top_block(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Top Block")

        ##################################################
        # Variables
        ##################################################
        cwd = os.getcwd()  # MODIFIED Gets current directory, allows any user name to specified when setting up UDOO
        log_path_username = cwd[:-4]  # MODIFIED
        self.gainRF = gainRF = 17.0
        self.gainIF = gainIF = 8.0
        self.prefix = prefix = log_path_username + "/CURRENT_DATA/"  # MODIFIED this varible from the standard top block
        self.gains = gains = "RF=" + str(gainRF) + " dB - " + "IF=" + str(gainIF) + " dB - "
        self.samp_rate = samp_rate = 2500000
        self.recfile = recfile = prefix + gains + datetime.now().strftime("%Y-%m-%d-T%H.%M.%S") + ' - ' + datetime.now().strftime('%s.%f')[:-3] + ".dat"
        self.freq = freq = 148.000000e6
        self.audio_rate = audio_rate = 48000 

        ##################################################
        # Blocks
        ##################################################
        self.timesync_pixstream_source_0 = timesync.pixstream_source(log_path_username + '/CURRENT_DATA')  # MODIFIED this variable too
        self.rational_resampler_xxx_0_0 = filter.rational_resampler_ccc(
                interpolation=audio_rate,
                decimation=2500000,
                taps=None,
                fractional_bw=None,
        )
        self.osmosdr_source_0 = osmosdr.source( args="numchan=" + str(1) + " " + 'airspy=0,pack=1' )
        self.osmosdr_source_0.set_sample_rate(samp_rate)
        self.osmosdr_source_0.set_center_freq(freq, 0)
        self.osmosdr_source_0.set_freq_corr(0, 0)
        self.osmosdr_source_0.set_dc_offset_mode(2, 0)
        self.osmosdr_source_0.set_iq_balance_mode(2, 0)
        self.osmosdr_source_0.set_gain_mode(False, 0)
        self.osmosdr_source_0.set_gain(gainRF, 0)
        self.osmosdr_source_0.set_if_gain(gainIF, 0)
        self.osmosdr_source_0.set_bb_gain(10, 0)
        self.osmosdr_source_0.set_antenna('', 0)
        self.osmosdr_source_0.set_bandwidth(0, 0)

        self.blocks_vector_sink_x_0_0 = blocks.vector_sink_f(7)
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, recfile, False)
        self.blocks_file_sink_0_0.set_unbuffered(False)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.osmosdr_source_0, 0), (self.rational_resampler_xxx_0_0, 0))
        self.connect((self.rational_resampler_xxx_0_0, 0), (self.blocks_file_sink_0_0, 0))
        self.connect((self.timesync_pixstream_source_0, 0), (self.blocks_vector_sink_x_0_0, 0))

    def get_gainRF(self):
        return self.gainRF

    def set_gainRF(self, gainRF):
        self.gainRF = gainRF
        self.osmosdr_source_0.set_gain(self.gainRF, 0)
        self.set_gains("RF=" + str(self.gainRF)+ " dB - " + "IF=" + str(self.gainIF) + " dB - ")

    def get_gainIF(self):
        return self.gainIF

    def set_gainIF(self, gainIF):
        self.gainIF = gainIF
        self.osmosdr_source_0.set_if_gain(self.gainIF, 0)
        self.set_gains("RF=" + str(self.gainRF)+ " dB - " + "IF=" + str(self.gainIF) + " dB - ")

    def get_prefix(self):
        return self.prefix

    def set_prefix(self, prefix):
        self.prefix = prefix
        self.set_recfile(self.prefix + self.gains + datetime.now().strftime("%Y-%m-%d-T%H.%M.%S") + ' - ' + datetime.now().strftime('%s.%f')[:-3] + ".dat")

    def get_gains(self):
        return self.gains

    def set_gains(self, gains):
        self.gains = gains
        self.set_recfile(self.prefix + self.gains + datetime.now().strftime("%Y-%m-%d-T%H.%M.%S") + ' - ' + datetime.now().strftime('%s.%f')[:-3] + ".dat")

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.osmosdr_source_0.set_sample_rate(self.samp_rate)

    def get_recfile(self):
        return self.recfile

    def set_recfile(self, recfile):
        self.recfile = recfile
        self.blocks_file_sink_0_0.open(self.recfile)

    def get_freq(self):
        return self.freq

    def set_freq(self, freq):
        self.freq = freq
        self.osmosdr_source_0.set_center_freq(self.freq, 0)

    def get_audio_rate(self):
        return self.audio_rate

    def set_audio_rate(self, audio_rate):
        self.audio_rate = audio_rate


def main(top_block_cls=top_block, options=None):

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
