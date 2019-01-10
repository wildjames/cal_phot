import numpy as np
import matplotlib.pyplot as plt

files = [
  'ASASSN-16kr_2018-10-12_r.calib',
  'ASASSN-16kr_2018-10-15_r.calib',
'ASASSN-16kr_2018-10-16_A_r.calib',
'ASASSN-16kr_2018-10-16_B_r.calib',
]
binsize = len(files)

print("Binning by {}".format(binsize))

master_ts = []
master_fl = []
master_fe = []

fig, ax = plt.subplots(2, figsize=[8, 6], sharex=True, sharey=True)
for file in files:
    data = np.loadtxt(file, delimiter=' ')
    ts = data[:,0]
    fl = data[:,1]
    fe = data[:,2]

    ax[1].step(ts, fl, label=file)

    master_ts += list(ts)
    master_fl += list(fl)
    master_fe += list(fe)


master_ts = np.array(master_ts)
master_fl = np.array(master_fl)
master_fe = np.array(master_fe)

sorted_args = np.argsort(master_ts)

master_ts = master_ts[sorted_args]
master_fl = master_fl[sorted_args]
master_fe = master_fe[sorted_args]


bin_ts = []
bin_fl = []
bin_fe = []

beg = 0
end = binsize
i = 0
while end < len(master_ts):
    beg = i*binsize
    i += 1
    end = i*binsize

    err = master_fe[beg:end] ** 2
    err = np.mean(err)
    err = np.sqrt(err)

    bin_ts.append(np.mean(master_ts[beg:end]))
    bin_fl.append(np.mean(master_fl[beg:end]))
    bin_fe.append(err)

ax[0].step(bin_ts, bin_fl)

plt.legend()
plt.tight_layout()
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

plt.close()