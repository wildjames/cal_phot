#!/usr/bin/env python3

import hipercam as hcam
import numpy as np
import matplotlib.pyplot as plt
import sys

xsplit = []

def on_click(event):
    global xsplit

    xsplit.append(event.xdata)

    print('Splitting at time {}'.format(xsplit[-1]))

try:
    lf = sys.argv[1]
    if lf[-4:] != '.log':
        lf += '.log'
except IndexError:
    print("Not enough arguments! Please give a logfile to split.")

data = hcam.hlog.Hlog.from_ascii(lf)

# Grab the header of that logfile too. Read lines until we run out of hashes, then skip the rest.
header = ''
datafile =  open(lf, 'r')

line = datafile.readline()
while line[0] == '#':
    header += line
    line = datafile.readline()

# grab the red red_lightcurve
red_lightcurve = data.tseries('1', '1') / data.tseries('1', '2')
xsplit = [red_lightcurve.t[0]]

fig, ax = plt.subplots()
fig.canvas.mpl_connect('button_press_event', on_click)
red_lightcurve.mplot(ax)
plt.show()

xsplit.append(red_lightcurve.t[-1])

if len(xsplit) == 2:
    print("No breakpoints entered! Nevermind.")
    exit()

# Now loop through the breakpoints
for i, x in enumerate(xsplit):
    if i > 0:
        xlo = xsplit[i-1]
        xhi = x
        print('{} - {}'.format(xsplit[i-1], x))
        slice = red_lightcurve[(red_lightcurve.t < xhi) & (red_lightcurve.t > xlo)]

        fig, ax = plt.subplots()
        slice.mplot(ax)
        plt.show()

        of = '{}_{}.log'.format(lf[:-4], i)
        print('Writing to {}...'.format(of))
        with open(of, 'w') as f:
            # [line] is currently a string containing the first line of data.
            print(line)
            f.write("# This logfile is the result of splitting the logfile {}\n# The following is the header from that file:\n".format(lf))
            f.write(header)
            f.write('# time, flx, fl_err, mask\n')

            print('xlo: {}  - xhi: {}'.format(xlo, xhi))
            # while the third column of the line is in range, write it to the file.
            while float(line.split(' ')[2]) <= xhi:
                f.write(line)
                line = datafile.readline()
                if line.strip() == '#':
                    f.write('#\n')
                    line = datafile.readline()
                elif line.strip() == '':
                    break
                    
datafile.close()