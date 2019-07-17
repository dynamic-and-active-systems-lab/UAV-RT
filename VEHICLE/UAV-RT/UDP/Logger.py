import sys
from datetime import datetime
import time
import os


class Logger(object):
    def __init__(self):
        global log_file_name
        self.terminal = sys.stdout
        cwd = os.getcwd()  # Gets current directory, allows any user name to specified at setup
        log_path_prefix = cwd + '/LOGS/'
        log_file_name = log_path_prefix + "udp_log - " + \
            datetime.now().strftime("%Y-%m-%d-T%H_%M_%S") + ".txt"
        self.log = open(log_file_name, "a", 0)

    def write(self, message):
        self.terminal.flush()
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


class Logger2(object):
    def __init__(self):
        global log_file_name
        self.terminal = sys.stderr
        self.log = open(log_file_name, "a", 0)

    def write(self, message):
        self.terminal.flush()
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


# class Logger3(object):
#     def __init__(self):
#         global log_file_name
#         self.terminal = sys.stdin
#         self.log = open(log_file_name, "a", 0)
#
#     def write(self, message):
#         self.terminal.write(message)
#         self.log.write(message)
#
#     def flush(self):
#         # this flush method is needed for python 3 compatibility.
#         # this handles the flush command by doing nothing.
#         # you might want to specify some extra behavior here.
#         pass
