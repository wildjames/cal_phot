import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        'nbin',
        help='The number of bins the data will be divided into',
        type=int,
        default=350
    )
    parser.add_argument(
        'phi',
        help='The phase range of the files will be truncated to within MIN MAX',
        type=float,
        nargs=2,
        default=[-np.inf, np.inf]
    )
    parser.add_argument(
        'files',
        help='The files to be binned together',
        type=str,
        nargs='+'
    )

    args = parser.parse_args()

    nbins = args.nbin
    phi_min, phi_max = args.phi
    files = args.files


    print("Binning into {} bins between phase {} and {}".format(nbins, phi_min, phi_max))

    # plotting area. Top axis will have output, bottom will have input.
    fig, ax = plt.subplots(2, figsize=[8, 6], sharex=True, sharey=True)
    ax[0].set_title('Binning down to {} points'.format(nbins))


    # Slap all the data in a dict
    master_data = {
        'ts': [],
        'fl': [],
        'fe': []
    }

    lab = True
    for f in files:
        print("-", f)
        data = np.loadtxt(f, delimiter=' ')

        master_data['ts'].extend(data[:,0])
        master_data['fl'].extend(data[:,1])
        master_data['fe'].extend(data[:,2])

        ax[1].step(data[:,0], data[:,1], label=f)
        if lab:
            ax[0].step(data[:,0], data[:,1], color='lightgrey',
            where='mid', label='Input Data')
            lab = False
        else:
            ax[0].step(data[:,0], data[:,1], where='mid', color='lightgrey')

    print("Input data have a {} points".format(len(master_data['ts'])))

    # Lets work with arrays, shall we?
    master_data = pd.DataFrame(master_data)

    # Sort the master DF
    master_data.sort_values('ts', ascending=True, inplace=True)

    # Enforce the phase limits
    master_data = master_data.loc[(master_data['ts'] >= phi_min) & (master_data['ts'] <= phi_max)]

    # I'll bin onto this axis...
    out_X = np.linspace(phi_min, phi_max, nbins)
    # I need to add an extra value to the top and tail of the linspace, or digitize gives me the values out to +/- inf
    sep = np.mean(out_X[1:] - out_X[:-1])
    out_X = np.insert(out_X, 0, out_X[0]-sep)
    out_X = np.append(out_X, out_X[-1]+sep)

    # Which bins do the values go into? The above defines bin edges!
    inds = np.digitize(master_data['ts'], out_X)
    master_data['bin'] = inds

    # These functions will be applied to the dataframe when taking means.
    def weighted_mean(x):
        return np.average(
            x,
            weights=master_data.loc[x.index, 'fe']**-2.
        )
    def calc_errs(x):
        d = master_data.loc[x.index, 'fe']
        return np.sqrt(np.sum(d**2)) / len(d)

    # This dict tells the agg what function to use when aggregating each column
    func_dict = {'fl': weighted_mean,
                'ts': weighted_mean,
                'fe': calc_errs}

    # Take a weighted mean of each bin.
    binned_master = master_data.groupby('bin').agg(weighted_mean)
    binned_master.sort_values('ts', ascending=True, inplace=True)

    # Enforce the phase limits.


    # Step plots are physically reasonable
    ax[0].step(
        binned_master['ts'], binned_master['fl'],
        color='black', where='mid', label='Binned lightcurve',
        zorder=10
    )
    ax[0].errorbar(
        binned_master['ts'], binned_master['fl'],
        binned_master['fe'],
        color='black', fmt='none', alpha=0.7,
        label=None,
        zorder=10
    )

    # Stretch the x-axis a little, to make it easier to see what's what
    extension = np.max([abs(0.1*phi_min), abs(0.1*phi_max)])
    lims = [phi_min - extension, phi_max + extension]
    ax[0].set_xlim(lims)

    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()
    # block=false so that the user can see the plot while deciding stuff
    plt.show(block=False)

    cont = input("Write to a file? y/n: ")
    if cont.lower().startswith('y'):
        oname = input("Enter a filename: ")
        if oname == '':
            exit()
        if not oname.endswith('.calib'):
            oname += '.calib'

        with open(oname, 'w') as f:
            f.write("# This file was produced by binning the following files down to {} points between phase {} and {}:\n".format(nbins, phi_min, phi_max))
            for cf in files:
                f.write("# {}\n".format(cf))
            f.write("#\n# phase, flux, error\n")
            for _, row in binned_master.iterrows():
                t = row['ts']
                fl = row['fl']
                fe = row['fe']
                f.write("{} {} {}\n".format(t, fl, fe))

        figname = oname.replace('.calib', '.pdf')
        plt.savefig(figname)

    plt.close()
