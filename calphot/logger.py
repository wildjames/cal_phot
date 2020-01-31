import sys
from io import TextIOWrapper

from ruamel.yaml import YAML

FNAME = "Calibration.txt"

class Logger(TextIOWrapper):
    '''
    Usage:
        sys.stdout = Logger()

    Takes all print outputs and mirrors them to a file.
    If the file, <inFile>, is supplied, this is copied into the head of the output.
    '''
    def __init__(self, logname, inFile=None):
        self.terminal = sys.stdout
        self.log = open(logname, "w")


        if inFile:
            self.log.write("#####################################    COPY OF INPUT FILE    #####################################\n")
            with open(inFile, 'r') as f:
                for line in f:
                    self.log.write(line)
            self.log.write("\n#####################################    END OF INPUT FILE     #####################################\n\n\n")
            self.log.write("\n##################################### BEGIN CALIBRATION OUTPUT #####################################\n")

    def write(self, message):
        self.terminal.write(message)
        if not 'Burning in' in message and not 'Sampling data' in message:
            self.log.write(message)

def header(inFile):
    fname = 'Calibration.txt'
    with open(fname, 'w') as o:
        o.write("#####################################    COPY OF INPUT FILE    #####################################\n")
        with open(inFile, 'r') as f:
            for line in f:
                o.write(line)
        o.write("\n#####################################    END OF INPUT FILE     #####################################\n\n\n")
        o.write("\n##################################### BEGIN CALIBRATION OUTPUT #####################################\n")


def printer(string, end='\n', fname=None, terminal=True):
    string = str(string)

    if fname == None:
        fname = FNAME

    with open(fname, 'a') as f:
        f.write(string+end)

    if terminal: print(string, end=end)
