import numpy as np
import matplotlib.pyplot as plt

def plot_all(files, oname):
    '''
    Takes a list of .calib files, and plots them all over each other as a step plot.

    Returns None
    '''
    if oname == '':
        print("  Saving overplotted eclipses to 'overplotted_eclipses.pdf'")
        oname = 'overplotted_eclipses'

    # filter the files so we only have green lightcurves
    files = [x for x in files if '_g.calib' in x]
    for f in files:
        print(f)

    
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                    '#f781bf', '#a65628', '#984ea3',
                    '#999999', '#e41a1c', '#dede00']

    fig, ax = plt.subplots()

    d = np.loadtxt(files[0], delimiter=' ')
    d = np.transpose(d)
    odata = d

    ax.step(d[0,:], d[1,:], label=files[0], color=CB_color_cycle[0])

    i = 1
    for f in files[1:]:
        d = np.loadtxt(f, delimiter=' ')
        d = np.transpose(d)
        ax.step(d[0,:], d[1,:], label=f, color=CB_color_cycle[i%9])
        # d[1,:] = d[1,:] / np.mean(d[1,:])
        odata = np.append(odata, d, axis=1)
        i+= 1

    plt.legend()
    plt.tight_layout()
    plt.savefig(oname+'.pdf')
    plt.show()

    return None