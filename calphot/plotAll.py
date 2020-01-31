from os import mkdir, path

import matplotlib.pyplot as plt
import numpy as np

from .logger import printer


def plot_all(files, oname, myLoc='.'):
    '''
    Takes a list of .calib files, and plots them all over each other as a step plot.

    Returns None
    '''
    plt.ion()
    printer("")
    if oname == '':
        printer("  Saving overplotted eclipses to 'overplotted_eclipses.pdf'")
        oname = 'overplotted_eclipses'

    oname = path.join(myLoc, 'MCMC_LIGHTCURVES', "FIGS", oname)
    directory = path.split(oname)[0]
    if not path.isdir(directory):
        mkdir(directory)

    print("Saving figures to {}".format(oname))


    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                      '#f781bf', '#a65628', '#984ea3',
                      '#999999', '#e41a1c', '#dede00']

    fig, ax = plt.subplots(figsize=[12,8])

    i = 0
    for f in files:
        d = np.loadtxt(f, delimiter=' ')
        d = np.transpose(d)
        ax.step(d[0,:], d[1,:], label=f, color=CB_color_cycle[i%9])
        i+= 1

    plt.title("Flux-calibrated eclispes")
    ax.set_ylabel('Flux, mJy')
    ax.set_xlabel("Phase")
    plt.legend()
    plt.tight_layout()
    plt.show()
    print("This figure will be saved as {}".format(oname+'.pdf'))
    input('Continue > ')
    print("---------------------------------------------------------")
    plt.savefig(oname+'.pdf')

    plt.close()
    plt.ioff()

    return None
