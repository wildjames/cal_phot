import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

files = [
    'SDSSJ0748_0_2017-02-24_r.calib',
    'SDSSJ0748_0_2018-02-04_r.calib',
    'SDSSJ0748_0_2018-12-17_r_1.calib',
    'SDSSJ0748_0_2018-12-17_r_2.calib',
]

nbins = 300
phi_min = -0.1
phi_max = +0.2

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
    data = np.loadtxt(f, delimiter=' ')

    master_data['ts'].extend(data[:,0])
    master_data['fl'].extend(data[:,1])
    master_data['fe'].extend(data[:,2])

    ax[1].step(data[:,0], data[:,1], label=f)
    if lab:
        ax[0].step(data[:,0], data[:,1], color='lightgrey', label='Input Data')
        lab = False
    else:
        ax[0].step(data[:,0], data[:,1], color='lightgrey')

# Lets work with arrays, shall we?
master_data = pd.DataFrame(master_data)

# Sort the master DF
master_data.sort_values('ts', ascending=True, inplace=True)

# I'll bin onto this axis...
out_X = np.linspace(phi_min, phi_max, nbins)
# I need to add an extra value to the top and tail of the linspace, or digitize gives me the values out to +/- inf
sep = np.mean(out_X[1:] - out_X[:-1])
out_X = np.insert(out_X, 0, out_X[0]-sep)
out_X = np.append(out_X, out_X[-1]+sep)

# Which bins do the values go into? The above defines bin edges!
inds = np.digitize(master_data['ts'], out_X)

# Initial dict
binned_master = pd.DataFrame(data=None,
    columns=master_data.columns)

# Need to go from 1 (the 0th index is junk, extending off to -inf), to the value of nbins (avoid OBOE)
for bin in range(1, nbins+1):
    # Which data?
    indeces = np.where(inds==bin)
    # Grab data
    sliced = master_data.iloc[indeces]
    # Mean data
    mean = sliced.mean()

    # The mean above does not handle errors correctly. Let's do that;
    # X = A+B; dX^2 = xA^2 + dB^2
    err = np.array(sliced['fe']) **2
    err = np.sum(err)
    err = np.sqrt(err)
    mean['fe'] = err

    # Add to binned df. Ignore index must be enabled for some reason.
    binned_master = binned_master.append(mean, ignore_index=True)

# Step plots are physically reasonable
ax[0].step(binned_master['ts'], binned_master['fl'], color='black', where='mid', label='Binned lightcurve')

# Don't print errorbars if they're too tiny?
all_mean = binned_master.mean()
if all_mean['fe'] > all_mean['fl']/10.:
    print("Signal/error ratio: {} -- Plotting error bars".format(all_mean['ff']/all_mean['fe']))
    ax[0].errorbar(binned_master['ts'], binned_master['fl'],
        yerr=binned_master['fe'], fmt='none', ecolor='black', capsize=0)

# Stretch the x-axis a little, to make it easier to see what's what
extension = np.max([abs(0.1*phi_min), abs(0.1*phi_max)])
lims = [phi_min - extension, phi_max + extension]
ax[0].set_xlim(lims)

ax[0].legend()
ax[1].legend()
plt.tight_layout()
plt.show()