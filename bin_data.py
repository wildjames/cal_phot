import numpy as np
import matplotlib.pyplot as plt
import os

files = [
# 'SDSSJ0748_0_2017-01-22_KG5.calib',

# 'SDSSJ0748_0_2017-02-15_KG5.calib',
# 'SDSSJ0748_0_2017-12-12_KG5.calib',

# 'SDSSJ0748_0_2017-02-14_g_1.calib',
# 'SDSSJ0748_0_2017-02-14_g_2.calib',

# 'SDSSJ0748_0_2018-12-16_g.calib',

# 'SDSSJ0748_0_2018-12-17_g.calib',

'SDSSJ0748_0_2017-02-24_r.calib',
'SDSSJ0748_0_2018-02-04_r.calib',
'SDSSJ0748_0_2018-12-17_r_1.calib',
'SDSSJ0748_0_2018-12-17_r_2.calib',

# 'SDSSJ0748_0_2018-02-05_r.calib',
]
binsize = 3
phi_min = -0.2
phi_max = +0.5


print("Binning by {}".format(binsize))
print("Trimming phase range to {:.2f}-{:.2f}".format(phi_min, phi_max))

# These will hold all the data
master_ts = []
master_fl = []
master_fe = []

# Get the data from the files
fig, ax = plt.subplots(2, figsize=[8, 6], sharex=True, sharey=True)
ax[0].set_title('Binning by {}'.format(binsize))
for file in files:
    data = np.loadtxt(file, delimiter=' ')
    ts = data[:,0]
    fl = data[:,1]
    fe = data[:,2]

    # Mask the data down to our window
    mask = (ts < phi_max) * (ts > phi_min)
    inrange = np.where(mask)
    ts = ts[inrange]
    fl = fl[inrange]
    fe = fe[inrange]

    print("The file {} has {} data".format(file, ts.shape))

    ax[0].step(ts, fl, alpha=0.3, color='black')
    ax[1].step(ts, fl, label=file)

    master_ts += list(ts)
    master_fl += list(fl)
    master_fe += list(fe)

# Lists are for wimps. Use arrays
master_ts = np.array(master_ts)
master_fl = np.array(master_fl)
master_fe = np.array(master_fe)

# Sort by time, so that when we bin adjacent data, they're actually relevant
sorted_args = np.argsort(master_ts)

master_ts = master_ts[sorted_args]
master_fl = master_fl[sorted_args]
master_fe = master_fe[sorted_args]

# Output lists
bin_ts = []
bin_fl = []
bin_fe = []

# Bin the data. This is probably not the most efficient way of doing this but hey ho
beg = 0
end = binsize
i = 0
while end < len(master_ts):
    beg = i*binsize
    i += 1
    end = i*binsize

    ## mean is added in quadrature
    err = master_fe[beg:end] ** 2
    err = np.mean(err)
    err = np.sqrt(err)

    bin_ts.append(np.mean(master_ts[beg:end]))
    bin_fl.append(np.mean(master_fl[beg:end]))
    bin_fe.append(err)

# Plot the bin result
ax[0].step(bin_ts, bin_fl)

plt.legend()
plt.tight_layout()
# block=false so that the user can see the plot while deciding stuff
plt.show(block=False)

print("The result has {} data points.".format(len(bin_ts)))
cont = input("Write to a file? y/n: ")
if cont.lower()[0] == 'y':
    oname = input("Enter a filename: ")
    oname += '.calib'
    with open(oname, 'w') as f:
        f.write("# This file was produced by binning the following files by {}:\n".format(binsize))
        for cf in files:
            f.write("# {}\n".format(cf))
        f.write("#\n# phase, flux, error\n")
        for t, fl, fe in zip(bin_ts, bin_fl, bin_fe):
            f.write("{} {} {}\n".format(t, fl, fe))

    figname = "../figs/"+oname.replace('.calib', '.pdf')
    # if os.path.isfile(figname):
    #     new_figname = figname.replace('.pdf', 'BIN{}.pdf'.format(binsize))
    #     print("Error! {} already exists. Saving as {} instead...".format(figname, new_figname))
    #     figname = new_figname
    plt.savefig(figname)

plt.close()
